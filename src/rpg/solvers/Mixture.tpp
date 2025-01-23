#ifndef RPG_MIXTURE_TPP
#define RPG_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <prdc/cuda/resources.h>

#include <cmath>

namespace Pscf { 
namespace Rpg
{

   // Constructor
   template <int D>
   Mixture<D>::Mixture()
    : ds_(-1.0),
      useBatchedFFT_(true),
      nParams_(0),
      meshPtr_(nullptr),
      hasStress_(false)
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
      read(in, "ds", ds_);
      readOptional(in, "useBatchedFFT", useBatchedFFT_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   /*
   * Create associations with a mesh, FFT, UnitCell, and WaveList object.
   */
   template <int D>
   void Mixture<D>::associate(Mesh<D> const & mesh, FFT<D> const & fft, 
                              UnitCell<D> const & cell, WaveList<D>& wavelist)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(fft.meshDimensions() == mesh.dimensions());
      UTIL_CHECK(cell.nParameter() > 0);

      // Assign internal pointers and variables
      meshPtr_ = &mesh;
      nParams_ = cell.nParameter();

      // Create associations for all blocks, set nParams in Polymer objects
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).setNParams(nParams_);
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).associate(mesh, fft, cell, wavelist);
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
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(mesh().size() > 0);

      // Allocate memory in all Block objects
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).allocate(ds_, useBatchedFFT_);
         }
      }

      // Allocate memory in all Solvent objects
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).allocate();
         }
      }
   }

   /*
   * Update solvers to account for new lattice parameters.
   */
   template <int D>
   void Mixture<D>::updateUnitCell()
   {
      for (int i = 0; i < nPolymer(); ++i) {
         for (int j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).updateUnitCell();
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
                            DArray< RField<D> > & cFields)
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
         polymer(i).compute(wFields);
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
         solvent(i).compute(wFields[monomerId]);

         // Add solvent contribution to relevant monomer concentration
         RField<D>& monomerField = cFields[monomerId];
         RField<D> const & solventField = solvent(i).concField();
         UTIL_CHECK(solventField.capacity() == nMesh);
         VecOp::addEqV(monomerField, solventField);
      }
      
      hasStress_ = false;
   }

   /*  
   * Compute total stress.
   */  
   template <int D>
   void Mixture<D>::computeStress()
   {   
      UTIL_CHECK(nParams_ > 0);

      int i, j;

      // Compute stress for each polymer.
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).computeStress();
      } 

      // Accumulate total stress 
      for (i = 0; i < nParams_; ++i) {
         stress_[i] = 0.0;
         for (j = 0; j < nPolymer(); ++j) {
            stress_[i] += polymer(j).stress(i);
         }   
      }
      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   /*
   * Is the ensemble canonical (i.e, closed for all species)?
   */
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

            blockCFields[sectionId] = solvent(i).concField();
         }
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif
