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
#include <pscf/chem/Monomer.h>
#include <pscf/chem/PolymerModel.h>
#include <util/containers/DArray.h>

namespace Pscf { 
namespace Rpg {

   using namespace Prdc;
   using namespace Prdc::Cuda;

   /*
   * Constructor
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

   /*
   * Destructor
   */
   template <int D>
   Mixture<D>::~Mixture()
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      // Read standard data for a mixture
      MixtureTmpl< PolymerT, SolventT >::readParameters(in);
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
   void Mixture<D>::associate(Mesh<D> const & mesh, 
                              FFT<D> const & fft, 
                              UnitCell<D> const & cell, 
                              WaveList<D> & waveList)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(fft.meshDimensions() == mesh.dimensions());
      UTIL_CHECK(cell.nParameter() > 0);

      // Assign member variables
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
      UTIL_CHECK(meshPtr_->size() > 0);
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

      // Allocate memory for all Solvent objects
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).allocate();
         }
      }

      clearUnitCellData();
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
            BlockT& block = polymer(i).block(j);
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
   void Mixture<D>::compute(DArray<FieldT> const & wFields, 
                            DArray<FieldT> & cFields,
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
      int i, j;

      // Clear all monomer concentration fields, check capacities
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == meshSize);
         UTIL_CHECK(wFields[i].capacity() == meshSize);
         VecOp::eqS(cFields[i], 0.0);
      }

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
            FieldT& monomerField = cFields[monomerId];
            FieldT const & blockField = polymer(i).block(j).cField();
            UTIL_CHECK(blockField.capacity() == meshSize);
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
         FieldT& monomerField = cFields[monomerId];
         FieldT const & solventField = solvent(i).cField();
         UTIL_CHECK(solventField.capacity() == meshSize);
         VecOp::addEqV(monomerField, solventField);
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

         // Compute stress for each polymer
         for (i = 0; i < nPolymer(); ++i) {
            polymer(i).computeStress();
         } 
   
         // Accumulate total stress 
         for (i = 0; i < nParam_; ++i) {
            for (j = 0; j < nPolymer(); ++j) {
               stress_[i] += polymer(j).stress(i);
            }   
         }

      }

      // Correct for partial occupation of the unit cell
      // Used in problems that contain a Mask, e.g., thin films
      for (i = 0; i < nParam_; ++i) {
         stress_[i] /= phiTot;
      }

      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   /*
   * Combine cFields for all blocks and solvents into one DArray
   */
   template <int D>
   void Mixture<D>::createBlockCRGrid(DArray<FieldT> & blockCFields) 
   const
   {
      int np = nSolvent() + nBlock();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(nMonomer() > 0);
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
         VecOp::eqS(blockCFields[i], 0.0);
      }

      #if 0
      UTIL_CHECK(blockCFields.capacity() == nBlock() + nSolvent());

      // Clear all monomer concentration fields, check capacities
      for (i = 0; i < np; ++i) {
         UTIL_CHECK(blockCFields[i].capacity() == nx);
         VecOp::eqS(blockCFields[i], 0.0);
      }
      #endif

      // Initialize section (block or solvent) index
      int sectionId = 0;

      // Process polymer species
      if (nPolymer() > 0) {
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               UTIL_CHECK(sectionId >= 0);
               UTIL_CHECK(sectionId < np);
               UTIL_CHECK(blockCFields[sectionId].capacity() == nx);
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

} // namespace Rpg
} // namespace Pscf
#endif
