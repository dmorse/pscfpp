#ifndef PRDC_CL_W_FIELDS_TPP
#define PRDC_CL_W_FIELDS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"
#include <prdc/fieldIo/cFieldIo.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <util/signal/Signal.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Cp {

   using namespace Util;
   using namespace Prdc;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D, class CFT, class FIT>
   WFields<D,CFT,FIT>::WFields()
    : fields_(),
      meshDimensions_(),
      meshSize_(0),
      nMonomer_(0),
      readUnitCellPtr_(nullptr),
      writeUnitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      signalPtr_(nullptr),
      isAllocated_(false),
      hasData_(false)
   {  signalPtr_ = new Signal<void>(); }

   /*
   * Destructor.
   */
   template <int D, class CFT, class FIT>
   WFields<D,CFT,FIT>::~WFields()
   {  delete signalPtr_; }

   /*
   * Create an association with a field Io (FIT) object.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::setFieldIo(FIT const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the unit cell that is modified by reading a field file.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::setReadUnitCell(UnitCell<D>& cell)
   {
      UTIL_CHECK(!readUnitCellPtr_);
      readUnitCellPtr_ = &cell;
   }

   /*
   * Set the unit cell that is written to a field file header.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::setWriteUnitCell(UnitCell<D> const & cell)
   {
      UTIL_CHECK(!writeUnitCellPtr_);
      writeUnitCellPtr_ = &cell;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::allocate(int nMonomer,
                                     IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(!hasData_);
      UTIL_CHECK(!isAllocated_);
      UTIL_CHECK(nMonomer > 0);

      // Store nMonomer and mesh dimensions
      nMonomer_ = nMonomer;
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }

      // Allocate arrays
      fields_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         fields_[i].allocate(meshDimensions);
      }

      isAllocated_ = true;
   }

   // Field Modification Functions

   /*
   * Set all field values.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::setFields(DArray<CFT> const & fields)
   {
      UTIL_CHECK(isAllocated_);

      // Update fields
      UTIL_CHECK(fields.capacity() == nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         UTIL_CHECK(fields[i].capacity() == meshSize_);
         assignField(fields_[i], fields[i]);
      }

      // Mark data as modified
      hasData_ = true;
      signal().notify();
   }

   // Read fields from a file

   /*
   * Write fields to an output stream.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::readFields(std::istream& in)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(readUnitCellPtr_);
      UTIL_CHECK(fieldIoPtr_);

      fieldIo().readFields(in, fields_, *readUnitCellPtr_);

      // Mark data as modified
      hasData_ = true;
      signal().notify();
   }

   /*
   * Read fields from a file, by filename.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::readFields(std::string const & filename)
   {
      std::ifstream file;
      fieldIo().fileMaster().openInputFile(filename, file);
      readFields(file);
      file.close();
   }

   // Write fields to file

   /*
   * Write fields to an output stream.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::writeFields(std::ostream& out) const
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(writeUnitCellPtr_);
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(hasData_);

      bool writeHeader = true;

      fieldIo().writeFields(out, fields_, *writeUnitCellPtr_, 
                            writeHeader);
   }

   /*
   * Write fields to a named file.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::writeFields(std::string const & filename) const
   {
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);
      writeFields(file);
      file.close();
   }

   // Signal accessor

   /*
   * Get the Signal<void> that is triggered by field modification.
   */
   template <int D, class CFT, class FIT>
   Signal<void>& WFields<D,CFT,FIT>::signal()
   {
      UTIL_CHECK(signalPtr_);
      return *signalPtr_;
   }

   // Private virtual function

   /*
   * Assignment operation for r-grid fields (CFT objects).
   *
   * Unimplemented virtual function - must be overridden by subclasses.
   */
   template <int D, class CFT, class FIT>
   void WFields<D,CFT,FIT>::assignField(CFT & lhs, CFT const & rhs) const
   {  UTIL_THROW("Unimplemented function WFields::assignRField"); }

} // namespace Cp
} // namespace Pscf
#endif
