#ifndef PSCF_HOST_D_ARRAY_H
#define PSCF_HOST_D_ARRAY_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace Pscf {

   // Forward declaration
   template <typename Data, typename T> class DeviceArray;
   class CUT;

   using namespace Util;

   /**
   * Template for dynamic array stored in host CPU memory.
   *
   * This class is provided as a convenience to allow the use of 
   * assigment (=) operators to copy data from device to host memory.
   * A HostArray<Data> stores data in a dynamically allocated array 
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
   class HostArray : public DArray<Data>
   {

   public:

      /**
      * Data type of each element.
      */
      using ValueType = Data;

      /**
      * Default constructor.
      */
      HostArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      HostArray(int capacity);

      HostArray(HostArray<Data> const & other) = default;

      /**
      * Copy constructor (copies from device to host).
      * 
      * \param other DeviceArray<Data,CUT> to be copied (input)
      */
      HostArray(DeviceArray<Data,CUT> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~HostArray();

      HostArray<Data>& operator = (HostArray<Data> const & other) = default;

      /**
      * Assignment operator, assign from a DeviceArray<Data,CUT>.
      *
      * Performs a deep copy from a RHS DeviceArray<Data,CUT> to this LHS
      * HostArray<D>, by copying the underlying data from device memory
      * to host memory.
      *
      * Preconditions: The RHS DeviceArray<Data,CUT> object must be allocated.  
      * If this LHS HostArray<D> is not allocated, the required memory
      * will be allocated before values are copied. Otherwise, if this LHS 
      * array is allocated on entry, capacites for LHS and RHS objects 
      * must be equal. 
      *
      * \throw Exception if the RHS array is not allocated on entry
      * \throw Exception if LHS and RHS have unequal nonzero capacities
      *
      * \param other DeviceArray<Data,CUT>  array on RHS of assignment (input)
      */
      HostArray<Data>& operator = (DeviceArray<Data,CUT> const & other);

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

   };

} // namespace Pscf

#include "DeviceArray.h"
#include "cudaErrorCheck.h"
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   /*
   * Default constructor.
   */
   template <typename Data>
   HostArray<Data>::HostArray()
    : DArray<Data>()
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   HostArray<Data>::HostArray(int capacity)
    : DArray<Data>()
   {  DArray<Data>::allocate(capacity); }

   /*
   * Copy constructor - deep copy DeviceArray from device to host.
   */
   template <typename Data>
   HostArray<Data>::HostArray(const DeviceArray<Data,CUT>& other)
    : DArray<Data>() 
   {  
      // Precondition - RHS array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data,CUT> must be allocated.");
      }

      DArray<Data>::allocate(other.capacity());
      cudaErrorCheck( cudaMemcpy(DArray<Data>::cArray(), other.cArray(), 
                                 DArray<Data>::capacity() * sizeof(Data), 
                                 cudaMemcpyDeviceToHost) );
   }

   /*
   * Destructor.
   */
   template <typename Data>
   HostArray<Data>::~HostArray()
   {} // DArray base class destructor will deallocate memory

   /*
   * Assignment from a DeviceArray<Data,CUT> RHS device array.
   */
   template <typename Data>
   HostArray<Data>& 
   HostArray<Data>::operator = (DeviceArray<Data,CUT> const & other)
   {
      // Precondition - RHS array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data,CUT> must be allocated.");
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
   void HostArray<Data>::copySlice(DeviceArray<Data,CUT> const & other,
                                    int beginId)
   {
      // Precondition - device array must be allocated
      if (!other.isAllocated()) {
         UTIL_THROW("RHS DeviceArray<Data,CUT> must be allocated.");
      }

      // Precondition - host array must be allocated
      if (!DArray<Data>::isAllocated()) {
         UTIL_THROW("LHS HostArray<Data> must be allocated.");
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
