#ifndef PSCF_HOST_D_ARRAY_CU_H
#define PSCF_HOST_D_ARRAY_CU_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>   // base class template
#include <pscf/backend/cuda/CUT.h>    // template argument

namespace Pscf {

   // Forward declarations
   template <typename Data, typename T> class HostArray;
   template <typename Data, typename T> class DeviceArray;

   using namespace Util;

   /**
   * Template for dynamic array stored in host CPU memory.
   *
   * This class template should be used in device-independent template
   * code in which an array is copied from host to device.  The
   * DeviceArray<Data,CPT> template defines an assignment operator
   * that assigns from a host array to a device array. 
   *
   * The "associate" member function allocates the array if not allocated
   * previously, or does nothing if the array is already allocated with
   * the correct capacity.
   *
   * The assignment (=) operator that copies a host array (RHS) to a
   * device array (LHS) performs a deep copy from host to device memory.
   * This function is a member function of the device array, defined by
   * the DeviceArray<Data,CUT> class template. 
   *
   * <b> Usage </b>: 
   *
   * Typical usage for backend-independent template code is shown 
   * below for host-to-device transfer to a long lived instance of 
   * DeviceArray<Data,CPT> named dArray from a shorter lived instance of 
   * HostArray<Data,CPT> named hArray. Here, the alias Data denotes the 
   * type of each array element.
   *
   * \code
   *    HostArray<Data,CPT> hArray;
   *    hArray.associate(dArray);
   *
   *    \\ ( Initialize data in hArray )
   *
   *    dArray = hArray
   *    hArray.dissociate();
   * \endcode
   *
   * Comments:
   *
   *   - In this specialization for a CUDA backend (T=CUT), the associate
   *     function allocates the host array, if not allocated previously,
   *     or does nothing if it is already allocated with the same 
   *     capacity as the device array. In the specialization for a C++
   *     backend (T=CPT), the associate function creates an association 
   *     that make the host array refer to memory owned by the device 
   *     array.
   *     
   *   - In this specialization for a CUDA backend (T=CUT), the 
   *     assignment (=) operator that assigns a RHS HostArray<Data,CPT> 
   *     to a LHS DeviceArray<Data,CPT> template copies all elements of
   *     an array from CPU host memory to GPU device memory. In the 
   *     corresponding specialization for a C++ backend (T=CPT), the
   *     assignment operator does nothing. 
   * 
   *   - In this specialization for a CUDA backend (T=CUT), the dissociate
   *     function does nothing. In the corresponding specialization for a
   *     C++ backend (T=CPT), this function destroys the association 
   *     between  the host and device arrays, by nullifying a pointer 
   *     held by the host array. 
   *
   *   - The host array may never be used to modify data after the 
   *     assignment operator and before the dissociate function is
   *     invoked. Doing so would modify data owned by the device array
   *     in CPU code (T=CPT) but would have no effect on data owned by
   *     the device array in GPU code (T=CUT), causing inconsistent
   *     behavior. To enforce this, it is good practice to invoke the 
   *     dissociate member function immediately after host-to-device 
   *     assignment, as shown above.
   *
   * \see Pscf::HostArray<Data,CPT>
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

      /**
      * Backend identifier class type.
      */
      using BackendIdClass = CUT;

      // Public member functions

      // Default constructor (default)
      HostArray() = default;

      /**
      * Allocating constructor.
      *
      * \param capacity  number of elements to allocate
      */
      HostArray(int capacity);

      // Copy constructor (deleted)
      HostArray(HostArray<Data,CUT> const & other) = delete;

      // Destructor (default).
      ~HostArray() = default;

      // Assignment (deleted).
      HostArray<Data,CUT>& 
      operator = (HostArray<Data,CUT> const & other) = delete;

      /**
      * Allocate if not allocated previously. 
      *
      * GPU specialization allocates host array with same dimensions as
      * the device array, unless this is already the case.
      *
      * \param deviceArray  device array (must be allocated on entry)
      */
      void associate(DeviceArray<Data,CUT>& deviceArray);

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
      * Release host array.
      *
      * GPU specialization (T=CUT) does nothing. The CPU specialization
      * (T=CPT) removes an association between the host array and some
      * other array. 
      */
      void dissociate()
      {}

      // Inherited public member functions (selected)
      using Array<Data>::operator [];
      using Array<Data>::cArray;
      using Array<Data>::isAllocated;
      using Array<Data>::capacity;
      using DArray<Data>::allocate;
      using DArray<Data>::deallocate;

   };

   // Explicit instantiation declarations
   extern template class HostArray<cudaReal,CUT>;
   extern template class HostArray<cudaComplex,CUT>;

} // namespace Pscf

#include <pscf/backend/cuda/DeviceArray.h>
#include <pscf/backend/cuda/cudaErrorCheck.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Allocating constructor.
   */
   template <typename Data>
   HostArray<Data,CUT>::HostArray(int capacity)
    : DArray<Data>(capacity) 
   {}

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
      UTIL_CHECK(isAllocated());
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
      UTIL_CHECK(beginId >= 0);
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
   * Allocate if not done previously. 
   */
   template <typename Data>
   void HostArray<Data,CUT>::associate(
                                 DeviceArray<Data,CUT>& deviceArray)
   {
      UTIL_CHECK(deviceArray.isAllocated());
      const int n = deviceArray.capacity();
      if (!isAllocated()) {
         allocate(deviceArray.capacity());
      }
      UTIL_CHECK(isAllocated());
      UTIL_CHECK(capacity() == n);
      // If this was allocated with capacity() == n, do nothing.
      // If this was allocated with capacity() != n, throw Exception.
   }

}
#endif
