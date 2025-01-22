#ifndef PSCF_HOST_D_ARRAY_H
#define PSCF_HOST_D_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   // Forward declaration of analogous container for data on the device.
   template <typename Data> class DeviceArray;

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
   * host memory.
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
      HostDArray<Data>& operator = (const DeviceArray<Data>& other);

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
   HostDArray<Data>::operator = (const DeviceArray<Data>& other)
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
      cudaMemcpy(DArray<Data>::cArray(), other.cArray(), 
                 DArray<Data>::capacity() * sizeof(Data), 
                 cudaMemcpyDeviceToHost);

      return *this;
   }

}
#endif
