#ifndef PSCF_HOST_D_ARRAY_CU_H
#define PSCF_HOST_D_ARRAY_CU_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>   // base class

namespace Pscf {

   // Forward declarations
   template <typename Data, typename T> class HostArray;
   template <typename Data, typename T> class DeviceArray;
   class CUT;

   using namespace Util;

   /**
   * Template for dynamic array stored in host CPU memory.
   *
   * This class is provided as a convenience to allow the use of 
   * assigment (=) operators to copy data from device to host memory.
   * A HostArray<Data,CUT> stores data in a dynamically allocated array 
   * in host CPU memory, whereas a DeviceArray<Data,CUT> stores analogous 
   * data in global GPU device memory. Each of these classes defines  
   * an assignment operation that allows assignment from the other, 
   * which silently copies the underlying arrays between device and 
   * host memory. Additionally, a method HostArray::copySlice is
   * provided, which populates a HostArray with a slice of a larger
   * DeviceArray.
   *
   * Otherwise, this class is identical to Util::DArray, with the
   * addition of an allocating constructor.
   *
   * \ingroup Pscf_Backend_Cuda_Module
   */
   template <typename Data>
   class HostArray<Data,CUT> : public DArray<Data>
   {

   public:

      /**
      * Data type of each element.
      */
      using ValueType = Data;

      // Default constructor (default)
      HostArray() = default;

      // Copy constructor (default)
      HostArray(HostArray<Data,CUT> const & other) = default;

      /**
      * Copy constructor (copies from device to host).
      * 
      * \param other DeviceArray<Data,CUT> to be copied (input)
      */
      HostArray(DeviceArray<Data,CUT> const & other);

      // Destructor (default).
      ~HostArray() = default;

      // Assignment (default)
      HostArray<Data,CUT>& 
      operator = (HostArray<Data,CUT> const & other) = default;

      /**
      * Assignment operator, assign from a DeviceArray<Data,CUT>.
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this LHS
      * HostArray<D>, by copying the underlying data from device memory
      * to host memory.
      *
      * Preconditions: The RHS DeviceArray<Data,CUT> object must be 
      * allocated.  If this LHS HostArray<D> is not allocated, the 
      * required memory will be allocated before values are copied. 
      * Otherwise, if this LHS array is allocated on entry, capacites 
      * for LHS and RHS objects must be equal. 
      *
      * \throw Exception if the RHS array is not allocated on entry
      * \throw Exception if LHS and RHS have unequal nonzero capacities
      *
      * \param other  DeviceArray<Data,CUT> on RHS of assignment (input)
      */
      HostArray<Data,CUT>& operator = (DeviceArray<Data,CUT> const & other);

      /**
      * Copy a slice of the data from a larger DeviceArray into this array.
      * 
      * This method will populate this HostArray with data from a slice
      * of a DeviceArray. The size of the slice is the capacity of this
      * HostArray (i.e., this entire array will be populated), and the
      * position of the slice within the DeviceArray is indicated by the
      * input parameter beginId. 
      * 
      * The capacity of the RHS DeviceArray must thus be greater than or
      * equal to the sum of beginId and the capacity of this HostArray.
      * 
      * \param other DeviceArray<Data,CUT>  object from which to copy slice
      * \param beginId  index of other array at which slice begins
      */
      void copySlice(DeviceArray<Data,CUT> const & other, int beginId);

      /**
      * Setup host array for use with a device array.
      *
      * GPU specialization allocates host array with same dimensions as
      * the device array, unless this is already the case.
      *
      * \param deviceArray  device array (must be allocated on entry)
      */
      void associate(DeviceArray<Data,CUT>& deviceArray);

      /**
      * Release host array.
      *
      * GPU specialization does nothing.
      */
      void dissociate()
      {}

      // Inherited public member functions (selected)
      using Array<Data>::capacity;
      using Array<Data>::isAllocated;
      using Array<Data>::operator [];
      using Array<Data>::cArray;
      using DArray<Data>::allocate;
      using DArray<Data>::deallocate;

   };

} // namespace Pscf

#include <pscf/backend/cuda/CUT.h> 
#include <pscf/backend/cuda/DeviceArray.h>
#include <pscf/backend/cuda/cudaErrorCheck.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Copy constructor - deep copy DeviceArray from device to host.
   */
   template <typename Data>
   HostArray<Data,CUT>::HostArray(DeviceArray<Data,CUT> const& other)
    : DArray<Data>() 
   {  
      // Precondition - RHS array must be allocated
      UTIL_CHECK(other.isAllocated());
      allocate(other.capacity());
      cudaErrorCheck( 
         cudaMemcpy(Array<Data>::data_, 
                    other.cArray(), 
                    capacity() * sizeof(Data), 
                    cudaMemcpyDeviceToHost) 
      );
   }

   /*
   * Assignment from a DeviceArray<Data,CUT> RHS device array.
   */
   template <typename Data>
   HostArray<Data,CUT>& 
   HostArray<Data,CUT>::operator = (DeviceArray<Data,CUT> const & other)
   {
      // Precondition - RHS array must be allocated
      UTIL_CHECK(other.isAllocated());

      // Allocate this if necessary 
      if (!isAllocated()) {
         allocate(other.capacity());
      } 

      // Require equal capacities
      UTIL_CHECK(capacity() == other.capacity());

      // Copy all elements
      cudaErrorCheck( 
         cudaMemcpy(Array<Data>::data_, 
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
   void HostArray<Data,CUT>::copySlice(DeviceArray<Data,CUT> const & other,
                                       int beginId)
   {
      // Preconditions 
      UTIL_CHECK(other.isAllocated());
      UTIL_CHECK(isAllocated());
      UTIL_CHECK(capacity() + beginId <= other.capacity());

      // Copy all elements
      cudaErrorCheck( 
         cudaMemcpy(Array<Data>::data_, 
                    other.cArray() + beginId, 
                    capacity() * sizeof(Data), 
                    cudaMemcpyDeviceToHost) 
      );
   }

   /*
   * Setup host array for use with device array. Allocate if necessary.
   */
   template <typename Data>
   void HostArray<Data,CUT>::associate(
                                 DeviceArray<Data,CUT>& deviceArray)
   {
      UTIL_CHECK(deviceArray.isAllocated());
      const int n = deviceArray.capacity();
      if (isAllocated() && capacity() != n) {
         deallocate();
      }
      if (!isAllocated()) {
         allocate(n);
      }
      // Note: If this was allocated with capacity() == n, nothing changes.
      UTIL_CHECK(capacity() == n);
   }

}
#endif
