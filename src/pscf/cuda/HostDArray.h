#ifndef PSCF_HOST_D_ARRAY_H
#define PSCF_HOST_D_ARRAY_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceArray.h"
#include "cudaErrorCheck.h"
#include <util/containers/DArray.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   /**
   * Template for dynamic array stored in host CPU memory.
   *
   * This class is provided as a convenience to allow the use of 
   * assigment (=) operators to copy data between corresponding 
   * containers that store array data in device vs. host memory. A 
   * HostDArray<Data> stores data in a dynamically allocated array in 
   * host CPU memory, whereas a DeviceArray<Data> stores analogous 
   * data in global GPU device memory. Each of these classes defines  
   * an assignment operation that allows assignment from the other, 
   * which silently copies the underlying arrays between device and 
   * host memory. Additionally, a method HostDArray::copySlice is
   * provided, which populates the HostDArray with a slice of the 
   * data from a larger DeviceArray.
   *
   * Otherwise, this class is identical to Util::DArray, with the
   * addition of an allocating constructor.
   *
   * \ingroup Pscf_Cuda_Module
   */
   template <typename Data>
   class HostDArray : public DArray<Data>
   {

   public:

      /**
      * Data type of each element.
      */
      typedef Data ElementType;

      /**
      * Default constructor.
      */
      HostDArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      HostDArray(int capacity);

      /**
      * Copy constructor. (Copies from any DArray or HostDArray).
      * 
      * \param other DArray<Data> to be copied (input)
      */
      HostDArray(DArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~HostDArray();

      /**
      * Assignment operator, assign from DeviceArray<Data> device array.
      *
      * Performs a deep copy from a RHS DeviceArray<Data> device array 
      * to this LHS HostDArray<D> host array, by copying the underlying 
      * C array from device memory to host memory.
      *
      * The RHS DeviceArray<Data> object must be allocated.  If this LHS 
      * HostDArray<D> is not allocated, the correct size block of memory 
      * will be allocated. Otherwise, if this LHS object is allocated on
      * entry, capacity values for LHS and RHS objects must be equal. 
      *
      * Exceptions are thrown if the RHS array is not allocated or if
      * both arrays are allocated but have unequal capacities.
      *
      * \param other DeviceArray<Data> array on RHS of assignment (input)
      */
      virtual 
      HostDArray<Data>& operator = (DeviceArray<Data> const & other);

      /**
      * Copy a slice of the data from a larger DeviceArray into this array.
      * 
      * This method will populate this HostDArray with data from a slice
      * of a DeviceArray. The size of the slice is the capacity of this
      * HostDArray (i.e., this entire array will be populated), and the
      * position of the slice within the DeviceArray is indicated by 
      * the input parameter beginId. 
      * 
      * Therefore, the capacity of the DeviceArray must be >= the 
      * capacity of this HostDArray plus beginId.
      * 
      * \param other DeviceArray<Data> object from which to copy slice
      * \param beginId  index of other array at which slice begins
      */
      void copySlice(DeviceArray<Data> const & other, int beginId);

   };

   /*
   * Default constructor.
   */
   template <typename Data>
   HostDArray<Data>::HostDArray()
    : DArray<Data>()
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   HostDArray<Data>::HostDArray(int capacity)
    : DArray<Data>()
   {  DArray<Data>::allocate(capacity); }

   /*
   * Copy constructor. (Copies from any DArray or HostDArray).
   */
   template <typename Data>
   HostDArray<Data>::HostDArray(const DArray<Data>& other)
    : DArray<Data>(other) // Use DArray base class copy constructor
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   HostDArray<Data>::~HostDArray()
   {} // DArray base class destructor will deallocate memory

   /*
   * Assignment from a DeviceArray<Data> RHS device array.
   */
   template <typename Data>
   HostDArray<Data>& 
   HostDArray<Data>::operator = (DeviceArray<Data> const & other)
   {
      // Precondition - RHS array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data> must be allocated.");
      }

      // Allocate this if necessary 
      if (!DArray<Data>::isAllocated()) {
         DArray<Data>::allocate(other.capacity());
      } 

      // Require equal capacities
      if (DArray<Data>::capacity() != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal capacity");
      }

      // Copy all elements
      cudaErrorCheck( cudaMemcpy(DArray<Data>::cArray(), other.cArray(), 
                                 DArray<Data>::capacity() * sizeof(Data), 
                                 cudaMemcpyDeviceToHost) );

      return *this;
   }

   /*
   * Copy a slice of the data from a larger DeviceArray into this array.
   */
   template <typename Data>
   void HostDArray<Data>::copySlice(DeviceArray<Data> const & other,
                                    int beginId)
   {
      // Precondition - device array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data> must be allocated.");
      }

      // Precondition - host array must be allocated
      if (!DArray<Data>::isAllocated()) {
         UTIL_THROW("LHS HostDArray<Data> must be allocated.");
      } 

      // Slice must not exceed the capacity of device array
      if (DArray<Data>::capacity() + beginId > other.capacity()) {
         UTIL_THROW("Slice must not exceed the capacity of device array.");
      }

      // Copy all elements
      cudaErrorCheck( cudaMemcpy(DArray<Data>::cArray(), 
                                 other.cArray() + beginId, 
                                 DArray<Data>::capacity() * sizeof(Data), 
                                 cudaMemcpyDeviceToHost) );
   }

}
#endif
