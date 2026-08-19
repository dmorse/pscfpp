#ifndef PSCF_DEVICE_MEMORY_H
#define PSCF_DEVICE_MEMORY_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/ReferenceCounter.h>

// Forward declaration
namespace Util { class CountedReference; }

namespace Pscf {

   using namespace Util;

   /**
   * Block of bare memory allocated on the GPU device.
   *
   * This class wraps a void* C array that is allocated in GPU device 
   * global memory, and also contains the number of bytes allocated.
   *
   * \ingroup Pscf_Cuda_Containers_Module
   */
   class DeviceMemory
   {

   public:

      /**
      * Default constructor.
      */
      DeviceMemory();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity size in bytes to allocate 
      */
      DeviceMemory(int capacity);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DeviceMemory();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if capacity parameter not positve
      * \throw Exception if the array is already allocated.
      *
      * \param capacity  number of byes to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array, if any.
      *
      * If no memory is allocated, do nothing. 
      */
      void deallocate();

      /**
      * Re-allocate if necessary to increase capacity.
      *
      * If the capacity parameter is greater than current capacity, allocate
      * or re-allocate to obtain a block of desired capacity. Otherwise, do
      * nothing. 
      *
      * When a previously allocated array is resized to increase capacity,
      * the old array is deallocated and all data is lost. 
      *
      * \throw Exception if capacity is not positve
      *
      * \param capacity  number of byes requird
      */
      void resize(int capacity);

      /**
      * Associate a reference with the reference counter.
      *
      * This should only be called within an associate function of a
      * container that creates an association with this object in 
      * which this object acts as the data owner. The CountedReference
      * should that is passed as a parameter should be owned by the
      * data user. On return, the ReferenceCounter that is owned by 
      * this DeviceMemory is incremented by one, and the CountedReference
      * contains a pointer to that ReferenceCounter.
      *
      * \param reference  reference to be included by reference counter
      */
      void addReference(CountedReference& reference);

      /**
      * Return pointer to underlying C array.
      */
      void* cArray() const;

      /**
      * Return allocated capacity.
      *
      * \return Number of elements allocated in array
      */
      int capacity() const;

      /**
      * Return true if the array has been allocated, false otherwise.
      */
      bool isAllocated() const;

   private:

      /// Pointer to a memory block on the GPU device.
      void* dataPtr_;

      /// Allocated size (capacity) of the array in bytes.
      int capacity_;

      /// Counter for references to data owned by this container.
      ReferenceCounter refCounter_;

   };

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   inline
   bool DeviceMemory::isAllocated() const
   {  return (bool)dataPtr_; }

   /*
   * Return allocated capacity.
   */
   inline 
   int DeviceMemory::capacity() const
   {  return capacity_; }

} // namespace Pscf
#endif
