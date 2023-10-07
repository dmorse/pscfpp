#ifndef PSPG_W_FIELD_CONTAINER_TPP
#define PSPG_W_FIELD_CONTAINER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFieldContainer.h"
#include <pspg/field/FieldIo.h>

namespace Pscf {
namespace Pspg
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
      nMonomer_(0),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false),
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
   * Set the stored value of nMonomer (this may only be called once).
   */
   template <int D>
   void WFieldContainer<D>::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
   }

   /*
   * Allocate memory for fields in r-grid format.
   */
   template <int D>
   void 
   WFieldContainer<D>::allocateRGrid(IntVec<D> const & meshDimensions)
   {
      UTIL_CHECK(nMonomer_ > 0);

      // If already allocated, deallocate.
      if (isAllocatedRGrid_) {
         deallocateRGrid();
      }
  
      // Store mesh dimensions
      meshDimensions_ = meshDimensions;
      meshSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshSize_ *= meshDimensions[i];
      }
  
      // Allocate arrays 
      rgrid_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].allocate(meshDimensions);
      }
      isAllocatedRGrid_ = true;

   }

   /*
   * De-allocate memory for fields in r-grid format
   */
   template <int D>
   void WFieldContainer<D>::deallocateRGrid()
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(nMonomer_ > 0);
      for (int i = 0; i < nMonomer_; ++i) {
         rgrid_[i].deallocate();
      }
      rgrid_.deallocate();
      meshDimensions_ = 0;
      meshSize_ = 0;
      isAllocatedRGrid_ = false;
   }
  
   /*
   * Allocate memory for fields in basis format.
   */
   template <int D>
   void WFieldContainer<D>::allocateBasis(int nBasis)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(nBasis > 0);

      // If already allocated, deallocate.
      if (isAllocatedBasis_) {
         deallocateBasis();
      }
 
      // Allocate 
      nBasis_ = nBasis;
      basis_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;

   }

   /*
   * De-allocate memory for fields in basis format.
   */
   template <int D>
   void WFieldContainer<D>::deallocateBasis()
   {
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(nMonomer_ > 0);
      for (int i = 0; i < nMonomer_; ++i) {
         basis_[i].deallocate();
      }
      basis_.deallocate();
      nBasis_ = 0;
      isAllocatedBasis_ = false;
   }
  
   /*
   * Allocate memory for all fields.
   */
   template <int D>
   void WFieldContainer<D>::allocate(int nMonomer, int nBasis, 
                                    IntVec<D> const & meshDimensions)
   {
      setNMonomer(nMonomer);
      allocateRGrid(meshDimensions);
      allocateBasis(nBasis);
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void 
   WFieldContainer<D>::setBasis(DArray< DArray<double> > const & fields)
   {
      UTIL_CHECK(isAllocatedBasis_);
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
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(fields.capacity() == nMonomer_);

      // Update system wFieldsRGrid
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize_, nBlocks, nThreads);
      for (int i = 0; i < nMonomer_; ++i) {
         //UTIL_CHECK(fields[i].capacity() == nBasis_)
         assignReal<<<nBlocks, nThreads>>>(rgrid_[i].cField(), 
                                           fields[i].cField(), 
                                           meshSize_);
         //rgrid_[i] = fields[i];
      }

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ =  isSymmetric;
   }

   /*
   * Set new w-field values, using unfoldeded array of r-grid fields.
   */
   template <int D>
   void WFieldContainer<D>::setRGrid(Field<cudaReal> & fields)
   {
      UTIL_CHECK(isAllocatedRGrid_);

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize_, nBlocks, nThreads);

      for (int i = 0; i < nMonomer_; i++) {
         assignReal<<<nBlocks, nThreads>>>(rgrid_[i].cField(), 
                                           fields.cField() + i*meshSize_, 
                                           meshSize_);
      }

      hasData_ = true;
      isSymmetric_ =  false;
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
      UTIL_CHECK(isAllocatedBasis());
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
      UTIL_CHECK(isAllocatedBasis());
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
      UTIL_CHECK(isAllocatedRGrid());
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
      UTIL_CHECK(isAllocatedRGrid());
      fieldIoPtr_->readFieldsRGrid(filename, rgrid_, unitCell);

      if (isSymmetric) {
         fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      }

      hasData_ = true;
      isSymmetric_ = isSymmetric;
   }

   /*
   * Set new w-field values, using array of r-grid fields as input.
   */
   template <int D>
   void WFieldContainer<D>::symmetrize()
   {
      UTIL_CHECK(hasData_);
      fieldIoPtr_->convertRGridToBasis(rgrid_, basis_);
      fieldIoPtr_->convertBasisToRGrid(basis_, rgrid_);
      isSymmetric_ = true;
   }

} // namespace Pspg
} // namespace Pscf
#endif
