#ifndef PRDC_CL_C_FIELDS_TPP
#define PRDC_CL_C_FIELDS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <prdc/fieldIo/cFieldIo.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
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
   CFields<D,CFT,FIT>::CFields()
    : fields_(),
      meshDimensions_(),
      meshSize_(0),
      nMonomer_(0),
      writeUnitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      isAllocated_(false),
      hasData_(false)
   {}

   /*
   * Destructor.
   */
   template <int D, class CFT, class FIT>
   CFields<D,CFT,FIT>::~CFields()
   {}

   /*
   * Create an association with a field Io (FIT) object.
   */
   template <int D, class CFT, class FIT>
   void CFields<D,CFT,FIT>::setFieldIo(FIT const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the unit cell that is written to a field file header.
   */
   template <int D, class CFT, class FIT>
   void CFields<D,CFT,FIT>::setWriteUnitCell(UnitCell<D> const & cell)
   {
      UTIL_CHECK(!writeUnitCellPtr_);
      writeUnitCellPtr_ = &cell;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D, class CFT, class FIT>
   void CFields<D,CFT,FIT>::allocate(int nMonomer,
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

   // Write fields to file

   /*
   * Write fields to an output stream.
   */
   template <int D, class CFT, class FIT>
   void CFields<D,CFT,FIT>::writeFields(std::ostream& out) const
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
   void CFields<D,CFT,FIT>::writeFields(std::string const & filename) const
   {
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);
      writeFields(file);
      file.close();
   }

   // Boolean flag setter

   /*
   * Set the hasData flag.
   */
   template <int D, class CFT, class FIT>
   void CFields<D,CFT,FIT>::setHasData(bool hasData) 
   {  hasData_ = hasData; }

} // namespace Cp
} // namespace Pscf
#endif
