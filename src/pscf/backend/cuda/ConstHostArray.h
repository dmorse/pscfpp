#ifndef PSCF_CONST_HOST_ARRAY_CU_H
#define PSCF_CONST_HOST_ARRAY_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstDArray.h>   // base class template
#include <pscf/backend/cuda/CUT.h>         // class template argument

namespace Pscf {

   // Forward declarations
   template <typename Data, typename T> class ConstHostArray;
   template <typename Data, typename T> class DeviceArray;

   using namespace Util;

   /**
   * Template for read-only dynamic array stored in host CPU memory.
   *
   * This class is derived from Util::ConstDArray<Data> and is almost
   * identical to this base class. The key difference is the addition of
   * specialized conversion constructor and assignment operators that 
   * allow construction and assignment from a DeviceArray<Data,CUT>.
   * Both of these functions copy an array from GPU device memory to 
   * host CPU memory.
   *
   * \ingroup Pscf_Backend_Cuda_Module
   */
   template <typename Data>
   class ConstHostArray<Data,CUT> : public ConstDArray<Data>
   {

   public:

      // Public types

      /**
      * Alias for data type of each array element.
      */
      using ValueType = Data;

      /**
      * Alias for backend identifier class.
      */
      using BackendIdClass = CUT;

      // Public member functions

      // Default constructor (default)
      ConstHostArray() = default;

      // Copy constructor (delete)
      ConstHostArray(ConstHostArray<Data,CUT> const & other) = delete;

      // Destructor (default)
      ~ConstHostArray() = default;

      // Assignment (delete)
      ConstHostArray<Data,CUT>&
      operator = (ConstHostArray<Data,CUT> const & other) = delete;

      /**
      * Conversion constructor (copies from device to host).
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this LHS
      * ConstHostArray<D>, by allocating this array with the required 
      * capacity and then copying the underlying array from device memory
      * to host memory.
      *
      * \throw Exception if the RHS array is not allocated on entry
      *
      * \param other DeviceArray<Data,CUT> to be copied (input)
      */
      ConstHostArray(DeviceArray<Data,CUT> const & other);

      /**
      * Assignment from a DeviceArray<Data,CUT>.
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this 
      * LHS ConstHostArray<D>, by first allocating if necessary, and then
      * copying the underlying data from device memory to host memory.
      *
      * \throw Exception if the RHS array is not allocated on entry
      * \throw Exception if LHS and RHS have unequal nonzero capacities
      *
      * \param other  DeviceArray<Data,CUT> on RHS of assignment (input)
      */
      ConstHostArray<Data,CUT>&
      operator = (DeviceArray<Data,CUT> const & other);

      /**
      * Copy a slice of data from a larger DeviceArray into this array.
      *
      * This function populate this ConstHostArray with data from a slice
      * of a DeviceArray. The size of the slice is equal to the capacity of
      * this LHS ConstHostArray, so that this entire array will be fully
      * populated. The position of the beginning of the slice within the 
      * RHS DeviceArray is indicated by input parameter beginId. The
      * capacity of the RHS DeviceArray must thus be greater than or equal
      * to the sum of beginId and the capacity of this LHS ConstHostArray.
      *
      * \throw Exception if RHS array is not allocated
      * \throw Exception if LHS array is not allocated
      * \throw Exception if beginId is negative
      * \throw Exception if the slice extends beyond the end of the RHS
      *
      * \param other DeviceArray<Data,CUT>  object from which to copy slice
      * \param beginId  index of other array at which slice begins
      */
      void copySlice(DeviceArray<Data,CUT> const & other, int beginId);

      /**
      * Release this array after finishing all reading.
      *
      * In template code that must work with either CPU or GPU, this should
      * be called after all code that requires read access to this array.
      * This GPU specialization (T=CUT) does nothing. The CPU specialization
      * (T=CPT) deletes the assocation with the device array.
      */
      void dissociate()
      {}

      // Inherited public member functions (selected)
      using ConstDArray<Data>::allocate;
      using ConstDArray<Data>::deallocate;
      using ConstDArray<Data>::operator =;
      using ConstArray<Data>::operator [];
      using ConstArray<Data>::cArray;
      using ConstArray<Data>::isAllocated;
      using ConstArray<Data>::capacity;
      using ConstArray<Data>::size;

   private:

      /**
      * Allocate if not allocated previously.
      *
      * If this is not allocated, allocate with a capacity equal to that
      * of the device array. 
      *
      * \throw Exception if allocated previously with wrong capacity
      *
      * \param deviceArray  device array (must be allocated on entry)
      */
      void associate(DeviceArray<Data,CUT> const & deviceArray);

   };

   // Explicit instantiation declarations
   extern template class ConstHostArray<cudaReal,CUT>;
   extern template class ConstHostArray<cudaComplex,CUT>;

} // namespace Pscf

#include "DeviceArray.h"
#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf 
{

   /*
   * Conversion constructor from DeviceArray<Data,CUT>.
   *
   * Allocate and perform a deep copy from device to host.
   */
   template <typename Data>
   ConstHostArray<Data,CUT>::ConstHostArray(
      DeviceArray<Data,CUT> const& other)
    : ConstDArray<Data>()
   {
      associate(other);

      // Copy data
      cudaErrorCheck( 
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray(),
                    capacity()*sizeof(Data),
                    cudaMemcpyDeviceToHost) 
      );
   }

   /*
   * Assignment from a DeviceArray<Data,CUT>.
   *
   * Allocate if necessary, and perform a deep copy from device to host.
   */
   template <typename Data>
   ConstHostArray<Data,CUT>&
   ConstHostArray<Data,CUT>::operator = (
                                  DeviceArray<Data,CUT> const & other)
   {
      // Allocate if not allocated previously, check capacities otherwise.
      associate(other);

      // Copy all elements
      cudaErrorCheck(
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray(),
                    capacity() * sizeof(Data),
                    cudaMemcpyDeviceToHost)
      );

      return *this;
   }

   /*
   * Copy a slice of the data from a larger DeviceArray into this array.
   */
   template <typename Data>
   void ConstHostArray<Data,CUT>::copySlice(
                                    DeviceArray<Data,CUT> const & other,
                                    int beginId)
   {
      // Preconditions
      UTIL_CHECK (other.isAllocated());
      UTIL_CHECK(isAllocated());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity() + beginId <= other.capacity());

      // Copy all elements
      cudaErrorCheck(
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray() + beginId,
                    capacity() * sizeof(Data),
                    cudaMemcpyDeviceToHost)
      );

   }

   // Private member function

   /*
   * Check allocation status, allocate if not allocated.
   */
   template <typename Data>
   void ConstHostArray<Data,CUT>::associate(
                               DeviceArray<Data,CUT> const & deviceArray)
   {
      UTIL_CHECK(deviceArray.isAllocated());
      const int n = deviceArray.capacity();
      if (!isAllocated()) {
         allocate(n);
      }
      UTIL_CHECK(capacity() == n);
      // If this was allocated with capacity() == n, do nothing.
      // If this was allocated with capacity() != n, throw Exception.
   }

}
#endif
