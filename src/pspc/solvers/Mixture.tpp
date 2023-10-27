#ifndef PSPC_MIXTURE_TPP
#define PSPC_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <prdc/cpu/RField.h>
#include <pscf/mesh/Mesh.h>

#include <cmath>

namespace Pscf {
namespace Pspc
{

   using namespace Prdc::Cpu;

   template <int D>
   Mixture<D>::Mixture()
    : ds_(-1.0),
      meshPtr_(0),
      unitCellPtr_(0),
      hasStress_(false)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl< Polymer<D>, Solvent<D> >::readParameters(in);
      read(in, "ds", ds_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   template <int D>
   void Mixture<D>::setDiscretization(Mesh<D> const & mesh,
                                      FFT<D> const & fft)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(mesh.dimensions() == fft.meshDimensions());

      // Save addresses of mesh
      meshPtr_ = &mesh;

      // Set discretization in space and s for all polymer blocks
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).setDiscretization(ds_, mesh, fft);
            }
         }
      }

      // Set spatial discretization for solvents
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).setDiscretization(mesh);
         }
      }

   }

   template <int D>
   void Mixture<D>::setupUnitCell(const UnitCell<D>& unitCell)
   {

      // Set association of unitCell to this Mixture
      unitCellPtr_ = &unitCell;

      // SetupUnitCell for all polymers
      if (nPolymer() > 0) {
         for (int i = 0; i < nPolymer(); ++i) {
            polymer(i).setupUnitCell(unitCell);
         }
      }

      hasStress_ = false;
   }


   /*
   * Reset statistical segment length for one monomer type.
   */
   template <int D>
   void Mixture<D>::setKuhn(int monomerId, double kuhn)
   {
      // Set new Kuhn length for relevant Monomer object
      monomer(monomerId).setKuhn(kuhn);

      // Update kuhn length for all blocks of this monomer type
      for (int i = 0; i < nPolymer(); ++i) {
         for (int j =  0; j < polymer(i).nBlock(); ++j) {
            if (monomerId == polymer(i).block(j).monomerId()) {
               polymer(i).block(j).setKuhn(kuhn);
            }
         }
      }
      hasStress_ = false;
   }

   /*
   * Compute concentrations (but not total free energy).
   */
   template <int D>
   void Mixture<D>::compute(DArray< RField<D> > const & wFields,
                            DArray< RField<D> > & cFields,
                            double phiTot)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(mesh().size() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      int meshSize = mesh().size();
      int nm = nMonomer();
      int i, j, k;

      // Clear all monomer concentration fields, check capacities
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == meshSize);
         UTIL_CHECK(wFields[i].capacity() == meshSize);
         for (j = 0; j < meshSize; ++j) {
            cFields[i][j] = 0.0;
         }
      }

      // Process polymer species
      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields, phiTot);
      }

      // Accumulate block contributions to monomer concentrations
      int monomerId;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            RField<D>& monomerField = cFields[monomerId];
            RField<D> const & blockField = polymer(i).block(j).cField();
            UTIL_CHECK(blockField.capacity() == meshSize);
            for (k = 0; k < meshSize; ++k) {
               monomerField[k] += blockField[k];
            }
         }
      }

      // Process solvent species
      // For each solvent, call compute and accumulate cFields
      for (i = 0; i < nSolvent(); ++i) {
         monomerId = solvent(i).monomerId();
         UTIL_CHECK(monomerId >= 0);
         UTIL_CHECK(monomerId < nm);

         // Compute solvent concentration
         solvent(i).compute(wFields[monomerId], phiTot);

         // Add solvent contribution to relevant monomer concentration
         RField<D>& monomerField = cFields[monomerId];
         RField<D> const & solventField = solvent(i).cField();
         UTIL_CHECK(solventField.capacity() == meshSize);
         for (k = 0; k < meshSize; ++k) {
            monomerField[k] += solventField[k];
         }

      }

      hasStress_ = false;
   }

   /*
   * Compute total stress for this mixture.
   */
   template <int D>
   void Mixture<D>::computeStress()
   {
      int i, j;

      // Initialize all stress components to zero
      for (i = 0; i < 6; ++i) {
         stress_[i] = 0.0;
      }

      if (nPolymer() > 0) {

         // Compute stress for all polymers, after solving MDE
         for (i = 0; i < nPolymer(); ++i) {
            polymer(i).computeStress();
         }
   
         // Accumulate stress for all the polymer chains
         for (i = 0; i < unitCellPtr_->nParameter(); ++i) {
            for (j = 0; j < nPolymer(); ++j) {
               stress_[i] += polymer(j).stress(i);
            }
         }

      }

      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   template <int D>
   bool Mixture<D>::isCanonical()
   {
      // Check ensemble of all polymers
      for (int i = 0; i < nPolymer(); ++i) {
         if (polymer(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Check ensemble of all solvents
      for (int i = 0; i < nSolvent(); ++i) {
         if (solvent(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Returns true if false was never returned
      return true;
   }

   /*
   * Combine cFields for each block (and solvent) into one DArray
   */
   template <int D>
   void
   Mixture<D>::createBlockCRGrid(DArray< RField<D> >& blockCFields) 
   const
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nBlock() + nSolvent() > 0);

      int np = nSolvent() + nBlock();
      int nx = mesh().size();
      int i, j;

      UTIL_CHECK(blockCFields.capacity() == nBlock() + nSolvent());

      // Clear all monomer concentration fields, check capacities
      for (i = 0; i < np; ++i) {
         UTIL_CHECK(blockCFields[i].capacity() == nx);
         for (j = 0; j < nx; ++j) {
            blockCFields[i][j] = 0.0;
         }
      }

      // Process polymer species
      int sectionId = -1;

      if (nPolymer() > 0) {

         // Write each block's r-grid data to blockCFields
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               sectionId++;

               UTIL_CHECK(sectionId >= 0);
               UTIL_CHECK(sectionId < np);
               UTIL_CHECK(blockCFields[sectionId].capacity() == nx);

               blockCFields[sectionId] = polymer(i).block(j).cField();
            }
         }
      }

      // Process solvent species
      if (nSolvent() > 0) {

         // Write each solvent's r-grid data to blockCFields
         for (i = 0; i < nSolvent(); ++i) {
            sectionId++;

            UTIL_CHECK(sectionId >= 0);
            UTIL_CHECK(sectionId < np);
            UTIL_CHECK(blockCFields[sectionId].capacity() == nx);

            blockCFields[sectionId] = solvent(i).cField();
         }
      }
   }

} // namespace Pspc
} // namespace Pscf
#endif
