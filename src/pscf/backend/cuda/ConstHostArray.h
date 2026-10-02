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
   * This class template should be used in backend-independent template
   * code in which data is copied from a device array to a host array.
   * This template is derived from ConstArray<Data>, and thus provides 
   * read-only access to the data.
   *
   * Conversion construction or assignment (operator =) from an instance 
   * of DeviceArray<Data,CUT> to a ConstHostArray<Data,CUT> allocates the
   * array if it is not allocated and then copies all elements of the 
   * array from global GPU device memory to host memory. Conversion
   * construction from a device array is equivalent to default 
   * construction followed by assignment.
   *
   * <b> Usage: </b>
   *
   * Typical usage is shown below for backend-indepent template code 
   * for device-to-host transfer from a longer lived instance of 
   * DeviceArray<Data,CUT> named dArray (the device array) to a 
   * shorter lived instance of ConstHostArray<Data,CUT> named hArray
   * (the host array). The type of each array element is denoted by Data.
   *
   * \code
   *    ConstHostArray<Data,CUT> hArray;
   *    hArray = dArray
   *
   *    // (Read and use the data in dArray)
   *
   *    hArray.dissociate();
   * \endcode
   * The default constructor and assignment operations may instead be 
   * combined into a single call of the conversion constructor, giving 
   * the shorter version
   * \code
   *    ConstHostArray<Data,CUT> hArray(dArray)
   *
   *    // Read and use the data in dArray
   *
   *    hArray.dissociate();
   * \endcode
   * The dissociate function does nothing in GPU code, but is used
   * for compatibility with corresponding CPU code, as discussed below.
   *
   * Comments:
   * 
   *  - In this specialization for a CUDA backend (T=CUT), the assignment
   *    operator and conversion constructor each perform a deep copy of 
   *    data from GPU device memory to CPU host memory. The corresponding 
   *    specialization for a GPU backend (T=CUT) would instead create a
   *    an association, i.e., a shallow copy, in which the host array
   *    owns a pointer to memory owned by the device array.
   *
   *  - In this specialization for a CUDA backend (T=CUT), the dissociate 
   *    function does nothing. It is used in backend-independent code to
   *    maintain compatibility with the syntax used by the specialization
   *    for a C++ backend (T=CPT), for which the dissociation member
   *    function releases the association with the device array, by 
   *    nullifying a pointer owned by the host array. 
   *
   *  - The dissociate function call must occur after all statements 
   *    that access elements of the host array. In the specialization 
   *    for a C++ backend (T=CPT), the host array has no access to data 
   *    after the dissociate function is called.
   *
   *  - If data that is copied from device to host needs to be modified
   *    on the host, it must first be copied to a container that provides
   *    read-write access. The ConstHostArray class template that is used
   *    for data transfer only provides read access.
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
