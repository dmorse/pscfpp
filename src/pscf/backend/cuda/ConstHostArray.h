#ifndef PSCF_CONST_HOST_ARRAY_CU_H
#define PSCF_CONST_HOST_ARRAY_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstDArray.h>   // base class

namespace Pscf {

   // Forward declarations
   template <typename Data, typename T> class ConstHostArray;
   template <typename Data, typename T> class DeviceArray;
   class CUT;

   using namespace Util;

   /**
   * Template for read-only dynamic array stored in host CPU memory.
   *
   * This class is derived from Util::ConstDArray<Data> and is almost
   * identical to this base class. The key difference is the addition of
   * specialized copy constructor and assignment operators that allow
   * construction and assignment from a DeviceArray<Data,CUT> object,
   * both of which copy data from GPU device memory to host CPU memory.
   *
   * \ingroup Pscf_Backend_Cuda_Module
   */
   template <typename Data>
   class ConstHostArray<Data,CUT> : public ConstDArray<Data>
   {

   public:

      /**
      * Data type of each element.
      */
      using ValueType = Data;

      // Default constructor (default)
      ConstHostArray() = default;

      // Copy constructor (deleted)
      ConstHostArray(ConstHostArray<Data,CUT> const & other) = delete;

      // Destructor
      ~ConstHostArray() = default;

      // Assignment (deleted)
      ConstHostArray<Data,CUT>&
      operator = (ConstHostArray<Data,CUT> const & other) = delete;

      /**
      * Copy constructor (copies from device to host).
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this LHS
      * ConstHostArray<D>, by copying the underlying data from device memory
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
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this LHS
      * ConstHostArray<D>, by copying the underlying data from device memory
      * to host memory.
      *
      * Preconditions: The RHS DeviceArray<Data,CUT> object must be
      * allocated.  If this LHS ConstHostArray<D> is not allocated, the
      * required memory will be allocated before values are copied.
      * Otherwise, if this LHS array is allocated on entry, capacites
      * for LHS and RHS objects must be equal.
      *
      * \throw Exception if the RHS array is not allocated on entry
      * \throw Exception if LHS and RHS have unequal nonzero capacities
      *
      * \param other  DeviceArray<Data,CUT> on RHS of assignment (input)
      */
      ConstHostArray<Data,CUT>&
      operator = (DeviceArray<Data,CUT> const & other);

      /**
      * Copy a slice of the data from a larger DeviceArray into this array.
      *
      * This method will populate this ConstHostArray with data from a slice
      * of a DeviceArray. The size of the slice is the capacity of this
      * ConstHostArray (i.e., this entire array will be populated), and the
      * position of the slice within the DeviceArray is indicated by the
      * input parameter beginId.
      *
      * The capacity of the RHS DeviceArray must thus be greater than or
      * equal to the sum of beginId and the capacity of this ConstHostArray.
      *
      * \param other DeviceArray<Data,CUT>  object from which to copy slice
      * \param beginId  index of other array at which slice begins
      */
      void copySlice(DeviceArray<Data,CUT> const & other, int beginId);

      /**
      * Setup host array for use with device array. Allocate if necessary.
      *
      * GPU specialization allocates host array with same dimensions as
      * the device array, unless this is already the case.
      *
      * \param deviceArray  device array (must be allocated on entry)
      */
      void associate(DeviceArray<Data,CUT> const & deviceArray);

      /**
      * Release host array.
      *
      * GPU specialization does nothing.
      */
      void dissociate()
      {}

      // Inherited public member functions (selected)
      using ConstArray<Data>::capacity;
      using ConstArray<Data>::size;
      using ConstArray<Data>::isAllocated;
      using ConstArray<Data>::cArray;
      using ConstArray<Data>::operator [];
      using ConstDArray<Data>::allocate;
      using ConstDArray<Data>::deallocate;
      using ConstDArray<Data>::operator =;

   };

} // namespace Pscf

#include "DeviceArray.h"
#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Copy constructor - deep copy DeviceArray from device to host.
   */
   template <typename Data>
   ConstHostArray<Data,CUT>::ConstHostArray(
      DeviceArray<Data,CUT> const& other)
    : ConstDArray<Data>()
   {
      // Precondition - RHS array must be allocated
      UTIL_CHECK(other.isAllocated());

      // If necessary, allocate this array
      if (!isAllocated()) {
          allocate(other.capacity());
      }

      // Copy data
      cudaErrorCheck( 
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray(),
                    capacity()*sizeof(Data),
                    cudaMemcpyDeviceToHost) 
      );
   }

   /*
   * Assignment from a DeviceArray<Data,CUT> RHS device array.
   */
   template <typename Data>
   ConstHostArray<Data,CUT>&
   ConstHostArray<Data,CUT>::operator = (
                                  DeviceArray<Data,CUT> const & other)
   {
      // Precondition - RHS array must be allocated
      UTIL_CHECK(other.isAllocated());

      // If necessary, allocate this array
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Require equal capacities
      if (capacity() != other.capacity()) {
         UTIL_THROW("Cannot assign arrays of unequal size");
      }

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
      UTIL_CHECK(capacity() + beginId <= other.capacity());

      // Copy all elements
      cudaErrorCheck(
         cudaMemcpy(ConstArray<Data>::data_,
                    other.cArray() + beginId,
                    capacity() * sizeof(Data),
                    cudaMemcpyDeviceToHost)
      );

   }

   /*
   * Setup host array for use with device array. Allocate if necessary.
   */
   template <typename Data>
   void ConstHostArray<Data,CUT>::associate(
                               DeviceArray<Data,CUT> const & deviceArray)
   {
      UTIL_CHECK(deviceArray.isAllocated());
      const int n = deviceArray.capacity();
      if (isAllocated() && capacity() != n) 
      { deallocate(); }
      if (!isAllocated()) {
         allocate(n);
      }
      // Note: If this was allocated with capacity() == n, nothing changes.
      UTIL_CHECK(capacity() == n);
   }

}
#endif
