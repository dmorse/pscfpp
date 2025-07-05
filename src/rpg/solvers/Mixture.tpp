#ifndef RPG_MIXTURE_TPP
#define RPG_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/chem/PolymerModel.h>

#include <util/containers/DArray.h>

//#include <cmath>

namespace Pscf { 
namespace Rpg
{

   /*
   * Constructor.
   */
   template <int D>
   Mixture<D>::Mixture()
    : stress_(),
      ds_(-1.0),
      meshPtr_(nullptr),
      nParam_(0),
      hasStress_(false),
      useBatchedFFT_(true)
   {  setClassName("Mixture"); }

   // Destructor
   template <int D>
   Mixture<D>::~Mixture()
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl< Polymer<D>, Solvent<D> >::readParameters(in);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);

      // Set ds_ parameter (only used in thread model)
      if (PolymerModel::isThread()) {
         read(in, "ds", ds_);
         UTIL_CHECK(ds_ > 0);
      } else {
         ds_ = 1.0;
      }

      // Optionally read useBatchedFFT boolean
      useBatchedFFT_ = true;
      readOptional(in, "useBatchedFFT", useBatchedFFT_);
   }

   /*
   * Create associations with a mesh, FFT, UnitCell, and WaveList object.
   */
   template <int D>
   void Mixture<D>::associate(Mesh<D> const & mesh, FFT<D> const & fft, 
                              UnitCell<D> const & cell, WaveList<D> & waveList)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(fft.meshDimensions() == mesh.dimensions());
      UTIL_CHECK(cell.nParameter() > 0);

      // Assign internal pointers and variables
      meshPtr_ = &mesh;
      nParam_ = cell.nParameter();

      // Create associations for all blocks, set nParams in Polymer objects
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            polymer(i).setNParams(nParam_);
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).associate(mesh, fft, cell, waveList);
            }
         }
      }

      // Create associations for all solvents
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).associate(mesh);
         }
      }

   }

   /*
   * Allocate internal data containers in all solvers. 
   */
   template <int D>
   void Mixture<D>::allocate()
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(mesh().size() > 0);
      UTIL_CHECK(ds_ > 0);

      // Allocate memory for all Block objects
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).allocate(ds_, useBatchedFFT_);
            }
         }
      }

      // Allocate memory in all Solvent objects
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).allocate();
         }
      }

      clearUnitCellData();
   }

   /*
   * Clear data that depends on lattice parameters in all solvers.
   */
   template <int D>
   void Mixture<D>::clearUnitCellData()
   {
      for (int i = 0; i < nPolymer(); ++i) {
         for (int j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).clearUnitCellData();
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

      int nMesh = mesh().size();
      int nm = nMonomer();
      int i, j;

      // Clear all monomer concentration fields
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == nMesh);
         UTIL_CHECK(wFields[i].capacity() == nMesh);
         VecOp::eqS(cFields[i], 0.0);
      }

      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields, phiTot);
      }

      // Accumulate monomer concentration fields
      int monomerId;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            RField<D>& monomerField = cFields[monomerId];
            RField<D>& blockField = polymer(i).block(j).cField();
            UTIL_CHECK(blockField.capacity() == nMesh);
            VecOp::addEqV(monomerField, blockField);
         }
      }
      
      // Process solvent species
      for (i = 0; i < nSolvent(); i++) {
         monomerId = solvent(i).monomerId();
         UTIL_CHECK(monomerId >= 0);
         UTIL_CHECK(monomerId < nm);

         // Compute solvent concentration
         solvent(i).compute(wFields[monomerId], phiTot);

         // Add solvent contribution to relevant monomer concentration
         RField<D>& monomerField = cFields[monomerId];
         RField<D> const & solventField = solvent(i).cField();
         UTIL_CHECK(solventField.capacity() == nMesh);
         VecOp::addEqV(monomerField, solventField);
      }
      
      hasStress_ = false;
   }

   /*  
   * Compute total stress.
   */  
   template <int D>
   void Mixture<D>::computeStress(double phiTot)
   {   
      UTIL_CHECK(nParam_ > 0);

      int i, j;

      // Compute stress for each polymer.
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).computeStress();
      } 

      // Accumulate total stress 
      for (i = 0; i < nParam_; ++i) {
         stress_[i] = 0.0;
         for (j = 0; j < nPolymer(); ++j) {
            stress_[i] += polymer(j).stress(i);
         }   
      }

      // Correct for partial occupation of the unit cell
      for (i = 0; i < nParam_; ++i) {
         stress_[i] /= phiTot;
      }

      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   /*
   * Is the ensemble canonical (i.e, closed for all species)?
   */
   template <int D>
   bool Mixture<D>::isCanonical() const
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
   void Mixture<D>::createBlockCRGrid(DArray< RField<D> > & blockCFields) 
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
         VecOp::eqS(blockCFields[i], 0.0);
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

} // namespace Rpg
} // namespace Pscf
#endif
