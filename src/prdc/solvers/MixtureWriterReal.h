#ifndef PRDC_MIXTURE_WRITER_REAL_H
#define PRDC_MIXTURE_WRITER_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}
#include <iostream>

namespace Pscf {
namespace Prdc {

   /**
   * A MixtureWriterReal writes fields stored by a Mixture to file.
   *
   * \ingroup Pscf_Prdc_Solvers_Module
   */
   template <int D, class MT, class FIT>
   class MixtureWriterReal 
   {

   public:

      // Public typename aliases

      using MixtureT = MT;
      using FieldIoT = FIT;
      using PolymerT = typename MT::PolymerT;
      using PropagatorT = typename MT::PropagatorT;
      using FieldT = typename MT::FieldT;

      /**
      * Constructor.
      */
      MixtureWriterReal();

      /**
      * Destructor.
      */
      ~MixtureWriterReal();

      /**
      * Create associations with related objects.
      */
      void associate(MixtureT const & mixture, 
                    FieldIoT const & fieldIo,
                    UnitCell<D> const & unitCell);

      /**
      * Set the isSymmetric flag true or false.
      *
      * This should be set true if the w fields used to solve the MDE for
      * all propagators is symmetric. A space group is written in field 
      * file headers only if a space group was declared in the parameter
      * file and isSymmetric is set true.
      */
      void setIsSymmetric(bool isSymmetric);

      /**
      * Write one slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write the complete propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name. The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

   protected:

      MixtureT const & mixture() const
      {  return *mixturePtr_; }

      FieldIoT const & fieldIo() const
      {  return *fieldIoPtr_; }

      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

   private:

      // Pointer to an associated MixtureT object 
      MixtureT const * mixturePtr_;

      // Pointer to an associated FieldIoT object 
      FieldIoT const * fieldIoPtr_;

      // Pointer to an associated UnitCell<D> object 
      UnitCell<D> const * unitCellPtr_;

      // Set true iff the w fields used in the MDE are symmetric
      bool isSymmetric_;

   };

} // namespace Prdc
} // namespace Pscf
#endif
