#ifndef RP_C_FIELDS_TPP
#define RP_C_FIELDS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <prdc/crystal/UnitCell.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   CFields<D,T>::CFields()
    : basis_(),
      rgrid_(),
      nMonomer_(0),
      writeUnitCellPtr_(nullptr),
      fieldIoPtr_(nullptr),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false),
      hasData_(false),
      isSymmetric_(false)
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D, class T>
   void
   CFields<D,T>::setFieldIo(FieldIo<D,T> const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Set the unit cell used for parameters written to a field header.
   */
   template <int D, class T>
   void CFields<D,T>::setWriteUnitCell(UnitCell<D> const & cell)
   {
      UTIL_CHECK(!writeUnitCellPtr_);
      writeUnitCellPtr_ = &cell;
   }

   /*
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D, class T>
   void CFields<D,T>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D, class T>
   void
   CFields<D,T>::allocateRGrid(IntVec<D> const & dimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!isAllocatedRGrid_);

      // Allocate arrays
      rgrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].allocate(dimensions);
      }
      isAllocatedRGrid_ = true;
   }

   /*
   * Allocate memory for fields in basis format.
   */
   template <int D, class T>
   void CFields<D,T>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(!isAllocatedBasis_);

      // Allocate
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * Allocate memory for all fields.
   */
   template <int D, class T>
   void
   CFields<D,T>::allocate(int nMonomer, int nBasis,
                                    IntVec<D> const & dimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(dimensions);
      allocateBasis(nBasis);
   }

   // Field output to file

   /*
   * Write fields to an output stream in basis format.
   */
   template <int D, class T>
   void CFields<D,T>::writeBasis(std::ostream& out) const
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(writeUnitCellPtr_);
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(hasData_);
      UTIL_CHECK(isSymmetric_);

      fieldIo().writeFieldsBasis(out, basis_, *writeUnitCellPtr_);
   }

   /*
   * Write fields to a file in basis format, by filename.
   */
   template <int D, class T>
   void CFields<D,T>::writeBasis(std::string filename) const
   {
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);
      writeBasis(file);
      file.close();
   }

   /*
   * Write fields to an output stream in real-space (r-grid) format.
   */
   template <int D, class T>
   void CFields<D,T>::writeRGrid(std::ostream& out) const
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(writeUnitCellPtr_);
      UTIL_CHECK(fieldIoPtr_);
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(hasData_);

      bool writeHeader = true;

      fieldIo().writeFieldsRGrid(out, rgrid_, *writeUnitCellPtr_, writeHeader,
                                 isSymmetric_);
   }

   /*
   * Write fields to a file in r-grid format, by filename.
   */
   template <int D, class T>
   void CFields<D,T>::writeRGrid(std::string filename) const
   {
      std::ofstream file;
      fieldIo().fileMaster().openOutputFile(filename, file);
      writeRGrid(file);
      file.close();
   }

   // Boolean flag setter functions

   // Set the hasData flag.
   template <int D, class T> inline
   void CFields<D,T>::setHasData(bool hasData)
   {
      hasData_ = hasData;
      if (!hasData_) {
         isSymmetric_ = false;
      }
   }

   // Set the isSymmetric flag.
   template <int D, class T> inline
   void CFields<D,T>::setIsSymmetric(bool isSymmetric)
   {
      UTIL_CHECK(hasData_);
      isSymmetric_ = isSymmetric;
   }

} // namespace Rp
} // namespace Pscf
#endif
