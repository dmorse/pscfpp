#ifndef PRDC_MIXTURE_REAL_TPP
#define PRDC_MIXTURE_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureReal.h"
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/chem/Monomer.h>
#include <pscf/chem/PolymerModel.h>
#include <util/containers/DArray.h>
#include <util/misc/ioUtil.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Prdc {

   /*
   * Constructor
   */
   template <int D, class PT, class ST>
   MixtureReal<D,PT,ST>::MixtureReal()
    : stress_(),
      ds_(-1.0),
      meshPtr_(nullptr),
      unitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      nParam_(0),
      hasStress_(false),
      isSymmetric_(false)
   {  setClassName("Mixture"); }

   /*
   * Destructor
   */
   template <int D, class PT, class ST>
   MixtureReal<D,PT,ST>::~MixtureReal()
   {}

   /*
   * Read all parameters and initialize.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::readParameters(std::istream& in)
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

   }

   /*
   * Create associations with Mesh, FFT, UnitCell, and WaveList objects.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::associate(Mesh<D> const & mesh,
                              FFTT const & fft,
                              UnitCell<D> const & cell,
                              WaveListT & waveList)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(mesh.size() > 0);
      UTIL_CHECK(fft.isSetup());
      UTIL_CHECK(mesh.dimensions() == fft.meshDimensions());
      UTIL_CHECK(cell.nParameter() > 0);

      // Assign member variables
      meshPtr_ = &mesh;
      unitCellPtr_ = &cell;
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

      // Create associations for all solvents (if any)
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).associate(mesh);
         }
      }

   }

   /*
   * Create an association with a FieldIo object, for file output.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::setFieldIo(FieldIoT const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate internal data containers in all solvers.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::allocate()
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(meshPtr_->size() > 0);
      UTIL_CHECK(ds_ > 0);

      // Allocate memory for all Block objects
      if (nPolymer() > 0) {
         allocateBlocks();
         #if 0
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).allocate(ds_);
            }
         }
         #endif
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
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::clearUnitCellData()
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
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::setKuhn(int monomerId, double kuhn)
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

   // Concentration field output

   /*
   * Compute concentrations (but not total free energy).
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::compute(DArray<FieldT> const & wFields,
                                      DArray<FieldT> & cFields,
                                      double phiTot)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(mesh().size() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      int nm = nMonomer();
      int monomerId;
      int i, j;

      // Clear all monomer concentration fields, check capacities
      for (i = 0; i < nm; ++i) {
         eqS(cFields[i], 0.0);
      }

      // Process polymer species
      if (nPolymer() > 0) {

         // Solve MDE for all polymers
         for (i = 0; i < nPolymer(); ++i) {
            polymer(i).compute(wFields, phiTot);
         }
   
         // Add block contributions to monomer concentrations
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               monomerId = polymer(i).block(j).monomerId();
               UTIL_CHECK(monomerId >= 0);
               UTIL_CHECK(monomerId < nm);
               FieldT& monomerField = cFields[monomerId];
               FieldT const & blockField = polymer(i).block(j).cField();
               addEqV(monomerField, blockField);
            }
         }

      }

      // Process solvent species
      if (nSolvent() > 0) {
      
         // Compute solvent concentrations
         for (i = 0; i < nSolvent(); ++i) {
            monomerId = solvent(i).monomerId();
            solvent(i).compute(wFields[monomerId], phiTot);
         }

         // Add solvent contributions to monomer concentrations
         for (i = 0; i < nSolvent(); ++i) {
            monomerId = solvent(i).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            FieldT& monomerField = cFields[monomerId];
            FieldT const & solventField = solvent(i).cField();
            addEqV(monomerField, solventField);
         }

      }

      isSymmetric_ = false;
      hasStress_ = false;
   }

   /*
   * Set the isSymmetric boolean variable true or false.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::setIsSymmetric(bool isSymmetric)
   {  isSymmetric_ = isSymmetric; }

   /*
   * Compute total SCFT stress for this mixture.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::computeStress(double phiTot)
   {
      UTIL_CHECK(nParam_ > 0);
      
      int i, j;

      // Initialize all stress components to zero
      stress_.clear();
      for (i = 0; i < nParam_; ++i) {
         stress_.append(0.0);
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

      // Correct for possible partial occupation of the unit cell.
      // Used in problems that contain a Mask, e.g., thin films.
      for (i = 0; i < nParam_; ++i) {
         stress_[i] /= phiTot;
      }

      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   // Concentration Field Output

   /*
   * Combine cFields for all blocks and solvents into one DArray
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::createBlockCRGrid(DArray<FieldT>& blockCFields)
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
         eqS(blockCFields[i], 0.0);
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

   /*
   * Output all concentration fields in real space (r-grid) format for
   * each block and solvent species to specified file.
   */
   template <int D, class PT, class ST>
   void 
   MixtureReal<D,PT,ST>::writeBlockCRGrid(std::string const & filename) 
   const
   {
      UTIL_CHECK(fieldIoPtr_);

      // Create and allocate array to hold c field data
      DArray<FieldT> blockCFields;
      int np = nSolvent() + nBlock();
      blockCFields.allocate(np);
      for (int i = 0; i < np; i++) {
         blockCFields[i].allocate(mesh().dimensions());
      }

      // Get c field data from the Mixture
      createBlockCRGrid(blockCFields);

      // Write block and solvent c field data to file
      fieldIo().writeFieldsRGrid(filename, blockCFields,
                                 unitCell(), isSymmetric_);
   }

   // Propagator field output

   /*
   * Write a specified slice of the propagator in r-grid format.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::writeQSlice(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId,
                                  int directionId, 
                                  int sliceId) const
   {
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < nPolymer());
      PolymerT const& polymerRef = polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymerRef.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      PropagatorT const &
         propagator = polymerRef.propagator(blockId, directionId);
      FieldT const & field = propagator.q(sliceId);
      fieldIo().writeFieldRGrid(filename, field, unitCell(), isSymmetric_);
   }

   /*
   * Write the last (tail) slice of the propagator in r-grid format.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::writeQTail(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId, 
                                  int directionId) const
   {
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < nPolymer());
      PolymerT const& polymerRef = polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymerRef.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      FieldT const &
         field = polymerRef.propagator(blockId, directionId).tail();
      fieldIo().writeFieldRGrid(filename, field, unitCell(), isSymmetric_);
   }

   /*
   * Write the entire propagator for a specified block and direction.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::writeQ(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId, 
                                  int directionId) const
   {
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < nPolymer());
      PolymerT const& polymerRef = polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymerRef.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      PropagatorT const&
           propagator = polymerRef.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);

      // Write header
      fieldIo().writeFieldHeader(file, 1, unitCell(), isSymmetric_);
      file << "mesh" << std::endl
           << "          " << mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
         file << "slice " << i << std::endl;
         fieldIo().writeFieldRGrid(file, propagator.q(i), unitCell(),
                                   hasHeader);
      }
      file.close();
   }

   /*
   * Write propagators for all blocks of all polymers to files.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::writeQAll(std::string const & basename)
   {
      UTIL_CHECK(fieldIoPtr_);
      std::string filename;
      int np, nb, ip, ib, id;
      np = nPolymer();
      for (ip = 0; ip < np; ++ip) {
         nb = polymer(ip).nBlock();
         for (ib = 0; ib < nb; ++ib) {
            for (id = 0; id < 2; ++id) {
               filename = basename;
               filename += "_";
               filename += toString(ip);
               filename += "_";
               filename += toString(ib);
               filename += "_";
               filename += toString(id);
               filename += ".rq";
               writeQ(filename, ip, ib, id);
            }
         }
      }
   }

   /*
   * Write stress to output stream.
   */
   template <int D, class PT, class ST>
   void MixtureReal<D,PT,ST>::writeStress(std::ostream& out) const
   {
      UTIL_CHECK(hasStress_);

      out << "stress:" << std::endl;

      for (int i = 0; i < nParam_; i++) {
         out << Int(i, 5)
            << "  "
            << Dbl(stress_[i], 18, 11)
            << std::endl;
      }
      out << std::endl;
   }

} // namespace Prdc
} // namespace Pscf
#endif
