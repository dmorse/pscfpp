#ifndef PSPC_W_FIELD_CONTAINER_TPP
#define PSPC_W_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.h"
#include <pspc/field/FieldIo.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   WFieldContainer<D>::WFieldContainer()
    : basis_(),
      rgrid_(),
      fieldIoPtr_(0),
      meshDimensions_(),
      meshSize_(0),
      nBasis_(0),
      isAllocated_(false),
      hasData_(false),
      isSymmetric_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   WFieldContainer<D>::~WFieldContainer()
   {}

   /*
   * Create an association with a FieldIo object.
   */
   template <int D>
   void WFieldContainer<D>::setFieldIo(FieldIo<D> const & fieldIo)
   {  fieldIoPtr_ = &fieldIo; }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void WFieldContainer<D>::allocate(int nMonomer, int nBasis, 
                                    IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(!isAllocated_);

      // Set mesh and basis dimensions
      nMonomer_ = nMonomer;
      nBasis_ = nBasis;
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate field arrays 
      basis_.allocate(nMonomer);
      rgrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         basis_[i].allocate(nBasis);
         rgrid_[i].allocate(meshDimensions);
      }

      isAllocated_ = true;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void WFieldContainer<D>::setBasis(DArray< DArray<double> > const & fields)
   {
      UTIL_CHECK(fields.capacity() == nMonomer_);

      // Update system wFields
      for (int i = 0; i < nMonomer_; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> &  w = basis_[i];
         UTIL_CHECK(f.capacity() == nBasis_);
         UTIL_CHECK(w.capacity() == nBasis_);
         for (int j = 0; j < nBasis_; ++j) {
            w[j] = f[j];
         }
      }

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Set new field values, using r-grid fields as inputs.
   */
   template <int D>
   void WFieldContainer<D>::setRGrid(DArray< RField<D> > const & fields,
                                    bool isSymmetric)
   {
      UTIL_CHECK(fields.capacity() == nMonomer_);

      // Update system wFieldsRGrid
      for (int i = 0; i < nMonomer_; ++i) {
         RField<D> const & f = fields[i];
         RField<D>& w = rgrid_[i];
         UTIL_CHECK(f.capacity() == meshSize_);
         UTIL_CHECK(w.capacity() == meshSize_);
         for (int j = 0; j < meshSize_; ++j) {
            w[j] = f[j];
         }
      }

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ =  isSymmetric;
   }

   /*
   * Read field component values from input stream, in symmetrized 
   * Fourier format.
   *
   * This function also computes and stores the corresponding
   * r-grid representation. On return, hasData and isSymmetric
   * are both true.
   */
   template <int D>
   void WFieldContainer<D>::readBasis(std::istream& in, 
                                      UnitCell<D>& unitCell)
   {
      UTIL_CHECK(isAllocated());
      fieldIoPtr_->readFieldsBasis(in, basis_, unitCell);

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Read field component values from file, in symmetrized 
   * Fourier format.
   *
   * This function also computes and stores the corresponding
   * r-grid representation. On return, hasData and isSymmetric
   * are both true.
   */
   template <int D>
   void WFieldContainer<D>::readBasis(std::string filename, 
                                      UnitCell<D>& unitCell)
   {
      UTIL_CHECK(isAllocated());
      fieldIoPtr_->readFieldsBasis(filename, basis_, unitCell);

      // Update system wFieldsRGrid
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);

      hasData_ = true;
      isSymmetric_ = true;
   }

   /*
   * Reads fields from an input stream in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the fields are known to be symmetric and so computes and stores
   * the corresponding basis components. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   * 
   * On return, hasData is true and the persistent isSymmetric flag 
   * defined by the class is set to the value of the isSymmetric 
   * input parameter.
   */
   template <int D>
   void WFieldContainer<D>::readRGrid(std::istream& in, 
                                      UnitCell<D>& unitCell, 
                                      bool isSymmetric)
   {
      UTIL_CHECK(isAllocated());
      fieldIoPtr_->readFieldsRGrid(in, rgrid_, unitCell);

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

   /*
   * Reads fields from a file in real-space (r-grid) format.
   *
   * If the isSymmetric parameter is true, this function assumes that 
   * the fields are known to be symmetric and so computes and stores
   * the corresponding basis components. If isSymmetric is false, it
   * only sets the values in the r-grid format.
   * 
   * On return, hasData is true and the persistent isSymmetric flag 
   * defined by the class is set to the value of the isSymmetric 
   * input parameter.
   */
   template <int D>
   void WFieldContainer<D>::readRGrid(std::string filename, 
                                      UnitCell<D>& unitCell, 
                                      bool isSymmetric)
   {
      UTIL_CHECK(isAllocated());
      fieldIoPtr_->readFieldsRGrid(filename, rgrid_, unitCell);

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

} // namespace Pspc
} // namespace Pscf
#endif
