#ifndef PSCF_HOST_D_ARRAY_TPP
#define PSCF_HOST_D_ARRAY_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostDArray.h"
#include "DeviceArray.h"

namespace Pscf {

   using namespace Util;

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
