#ifndef RPC_MIXTURE_TPP
#define RPC_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <prdc/cpu/RField.h>

#include <prdc/cpu/FFT.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/chem/Monomer.h>
#include <pscf/chem/PolymerModel.h>

#include <util/containers/DArray.h>

#include <cmath>

namespace Pscf {
namespace Rpc
{

   using namespace Prdc::Cpu;

   template <int D>
   Mixture<D>::Mixture()
    : stress_(),
      ds_(-1.0),
      meshPtr_(nullptr),
      nParam_(0),
      hasStress_(false)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      // Read majority of mixture
      MixtureTmpl< Polymer<D>, Solvent<D> >::readParameters(in);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);

      // Read ds parameter
      if (PolymerModel::isThread()) {
         read(in, "ds", ds_);
         UTIL_CHECK(ds_ > 0);
      } else {
         ds_ = 1.0;
      }

   }

   /*
   * Create associations with mesh, fft, and unit cell.
   */
   template <int D>
   void Mixture<D>::associate(Mesh<D> const & mesh,
                              FFT<D> const & fft,
                              UnitCell<D> const & cell)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(mesh.dimensions() == fft.meshDimensions());
      UTIL_CHECK(cell.nParameter() > 0);

      // Save addresses of mesh and unit cell
      meshPtr_ = &mesh;
      nParam_ = cell.nParameter();

      // Create associations for all polymer blocks
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            polymer(i).setNParams(nParam_);
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).associate(mesh, fft, cell);
            }
         }
      }

      // Create associations for all solvents (if any)
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).associate(mesh);
         }
      }

   }


   /*
   * Allocate internal data containers for all solvers.
   */
   template <int D>
   void Mixture<D>::allocate()
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(meshPtr_->size() > 0);
      UTIL_CHECK(ds_ > 0);

      // Allocate memory for all Block objects
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).allocate(ds_);
            }
         }
      }

      // Set spatial discretization for solvents
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).allocate();
         }
      }

   }

   /*
   * Clear all data that depends on the unit cell parameters.
   */
   template <int D>
   void Mixture<D>::clearUnitCellData()
   {
      if (nPolymer() > 0) {
         for (int i = 0; i < nPolymer(); ++i) {
            polymer(i).clearUnitCellData();
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
            Block<D>& block = polymer(i).block(j);
            if (monomerId == block.monomerId()) {
               block.setKuhn(kuhn);
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
   void Mixture<D>::computeStress(double phiTot)
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
         for (i = 0; i < nParam_; ++i) {
            for (j = 0; j < nPolymer(); ++j) {
               stress_[i] += polymer(j).stress(i);
            }
         }

      }

      // Correct for possible partial occupation of the unit cell.
      // Used in problems that contain a Mask, e.g., thin films.
      for (i = 0; i < nParam_; ++i) {
         stress_[i] /= phiTot;
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
   * Combine cFields for all blocks and solvents into one DArray
   */
   template <int D>
   void
   Mixture<D>::createBlockCRGrid(DArray< RField<D> >& blockCFields) 
   const
   {
      int np = nSolvent() + nBlock();
      UTIL_CHECK(np > 0);
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      int i, j;

      // Check allocation of blockCFields, allocate if necessary
      // Initialize all concentration values to zero
      if (!blockCFields.isAllocated()) {
         blockCFields.allocate(np);
      }
      UTIL_CHECK(blockCFields.capacity() == np);
      for (i = 0; i < np; ++i) {
         if (!blockCFields[i].isAllocated()) {
            blockCFields[i].allocate(mesh().dimensions());
         }
         UTIL_CHECK(blockCFields[i].capacity() == nx);
         for (j = 0; j < nx; ++j) {
            blockCFields[i][j] = 0.0;
         }
      }

      // Initialize section (block or solvent) index
      int sectionId = 0;

      // Process polymer species
      if (nPolymer() > 0) {
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               UTIL_CHECK(sectionId >= 0);
               UTIL_CHECK(sectionId < np);
               // Copy block r-grid c-field to blockCFields
               blockCFields[sectionId] = polymer(i).block(j).cField();
               sectionId++;
            }
         }
      }
      UTIL_CHECK(sectionId == nBlock());

      // Process solvent species
      if (nSolvent() > 0) {
         for (i = 0; i < nSolvent(); ++i) {
            UTIL_CHECK(sectionId >= 0);
            UTIL_CHECK(sectionId < np);
            // Copy solvent r-grid c-field to blockCFields
            blockCFields[sectionId] = solvent(i).cField();
            sectionId++;
         }
      }
      UTIL_CHECK(sectionId == np);

   }

} // namespace Rpc
} // namespace Pscf
#endif
