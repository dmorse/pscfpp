#ifndef PRDC_MIXTURE_WRITER_REAL_TPP 
#define PRDC_MIXTURE_WRITER_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureWriterReal.h"
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <util/misc/ioUtil.h>
#include <util/misc/FileMaster.h>

#include <string>
#include <unistd.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D, class MT, class FIT>
   MixtureWriterReal<D,MT,FIT>::MixtureWriterReal()
    : mixturePtr_(nullptr),
      fieldIoPtr_(nullptr),
      unitCellPtr_(nullptr),
      isSymmetric_(false)
   {}

   /*
   * Destructor.
   */
   template <int D, class MT, class FIT>
   MixtureWriterReal<D,MT,FIT>::~MixtureWriterReal()
   {}

   /*
   * Create associations with a Mixture, FieldIo and a UnitCell 
   */
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::associate(
                                  MixtureT const & mixture,
                                  FieldIoT const & fieldIo,
                                  UnitCell<D> const & unitCell)
   {
      mixturePtr_ = &mixture;
      fieldIoPtr_ = &fieldIo;
      unitCellPtr_ = &unitCell;
   }

   /*
   * Set the isSymmetric boolean variable true or false.
   */
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::setIsSymmetric(bool isSymmetric)
   {  isSymmetric_ = isSymmetric; }

   /*
   * Write a specified slice of the propagator in r-grid format.
   */
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::writeQSlice(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId,
                                  int directionId, 
                                  int sliceId) const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture().nPolymer());
      PolymerT const& polymer = mixture().polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      PropagatorT const &
         propagator = polymer.propagator(blockId, directionId);
      FieldT const & field = propagator.q(sliceId);
      fieldIo().writeFieldRGrid(filename, field, unitCell(), isSymmetric_);
   }

   /*
   * Write the last (tail) slice of the propagator in r-grid format.
   */
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::writeQTail(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId, 
                                  int directionId) const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture().nPolymer());
      PolymerT const& polymer = mixture().polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      FieldT const &
         field = polymer.propagator(blockId, directionId).tail();
      fieldIo().writeFieldRGrid(filename, field, unitCell(), isSymmetric_);
   }

   /*
   * Write the entire propagator for a specified block and direction.
   */
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::writeQ(
                                  std::string const & filename,
                                  int polymerId, 
                                  int blockId, 
                                  int directionId) const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture().nPolymer());
      PolymerT const& polymer = mixture().polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      UTIL_CHECK(unitCell().isInitialized());
      PropagatorT const&
           propagator = polymer.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);

      // Write header
      fieldIo().writeFieldHeader(file, 1, unitCell(), isSymmetric_);
      file << "mesh" << std::endl
           << "          " << fieldIo().mesh().dimensions() << std::endl
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
   template <int D, class MT, class FIT>
   void MixtureWriterReal<D,MT,FIT>::writeQAll(std::string const & basename)
   {
      std::string filename;
      int np, nb, ip, ib, id;
      np = mixture().nPolymer();
      for (ip = 0; ip < np; ++ip) {
         nb = mixture().polymer(ip).nBlock();
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

} // namespace Prdc
} // namespace Pscf
#endif
