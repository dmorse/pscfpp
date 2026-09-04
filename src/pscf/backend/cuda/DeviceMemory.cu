/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceMemory.h"
#include "cudaErrorCheck.h"
#include <util/misc/CountedReference.h>
#include <util/global.h>
#include <cuda_runtime.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   DeviceMemory::DeviceMemory()
    : dataPtr_(nullptr),
      capacity_(0),
      refCounter_()
   {}

   /*
   * Allocating constructor.
   */
   DeviceMemory::DeviceMemory(int capacity)
    : dataPtr_(nullptr),
      capacity_(0),
      refCounter_()
   {  allocate(capacity); }

   /*
   * Destructor.
   */
   DeviceMemory::~DeviceMemory()
   {
      if (isAllocated()) {
         if (refCounter_.hasRefs()) {
            std::cout << "Error - destruction of DeviceMemory with"
                      << "external references" << std::endl;
         }
         cudaFree(dataPtr_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   void DeviceMemory::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a DeviceMemory");
      }
      cudaErrorCheck( cudaMalloc( &dataPtr_, capacity ) );
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   */
   void DeviceMemory::deallocate()
   {
      UTIL_CHECK(!refCounter_.hasRefs());
      if (isAllocated()) {
         cudaErrorCheck( cudaFree(dataPtr_) );
      }
      dataPtr_ = nullptr;
      capacity_ = 0;
   }

   /*
   * Reallocate if necessary to increase capacity.
   */
   void DeviceMemory::resize(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to resize DeviceMemory with capacity <= 0");
      }
      if (capacity > capacity_) {
         deallocate();
         allocate(capacity);
      }
   }

   /*
   * Associate a CountedReference with the reference counter member.
   */
   void DeviceMemory::addReference(CountedReference& ref)
   {  ref.associate(refCounter_); } 

   /*
   * Get a pointer to the underlying C array.
   */
   void* DeviceMemory::cArray() const
   {
      UTIL_CHECK(dataPtr_);
      return dataPtr_;
   }

} // namespace Pscf
