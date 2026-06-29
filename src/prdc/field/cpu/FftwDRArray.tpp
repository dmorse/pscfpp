#ifndef PRDC_CPU_FFTW_DR_ARRAY_TPP
#define PRDC_CPU_FFTW_DR_ARRAY_TPP

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/cpu/FftwDRArray.h>
#include <util/misc/Memory.h> 

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {


   // Non-inline member function definitions

   /*
   * Default constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray()
    : Array<Data>()
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray(int capacity)
    : Array<Data>()
   {  allocate(capacity); }

   /*
   * Copy constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray(FftwDRArray<Data> const & other)
    : Array<Data>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other FftwDRArray must be allocated.");
      }
      allocate(other.capacity_);
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Destructor.
   */
   template <typename Data>
   FftwDRArray<Data>::~FftwDRArray()
   {
      if (data_) {
         if (ref_.isAssociated()) {
            ref_.dissociate();
         } else {
            if (refCounter_.hasRefs()) {
               int nRef = refCounter_.nRef();
               std::cout
                  << std::endl
                  << "Error: Destroying FftwDRArray that is referenced by "
                  << nRef << " other(s)" << std::endl;
            }
            try {
               fftw_free(data_);
               Memory::sub<Data>(capacity_);
            } catch (...) {
               std::cout
                  << std::endl
                  << "Error in deallocation in FftwDRArray destructor";
            }
         }
      }
      data_ = nullptr;
      capacity_ = 0;
   }

   /*
   * Assignment from another FftwDRArray<Data> (deep copy).
   */
   template <typename Data>
   FftwDRArray<Data>&
   FftwDRArray<Data>::operator = (FftwDRArray<Data> const & other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - other array must be allocated
      UTIL_CHECK (other.isAllocated());

      // If this array is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Copy elements
      UTIL_CHECK(capacity_ == other.capacity_);
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Assignment from an Array<Data> (deep copy).
   */
   template <typename Data>
   FftwDRArray<Data>&
   FftwDRArray<Data>::operator = (Array<Data> const & other)
   {
      // Check for self assignment
      if (dynamic_cast< Array<Data>* >(this) == &other) return *this;

      // Precondition - other array must be allocated
      UTIL_CHECK(other.capacity() > 0);

      // If this array is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Copy all elements
      UTIL_CHECK(capacity_ == other.capacity());
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate an underlying C array, which this container then owns.
   */
   template <typename Data>
   void FftwDRArray<Data>::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to allocate a FftwDRArray that has data.");
      }
      data_ = (Data*) fftw_malloc(sizeof(Data) * capacity);
      capacity_ = capacity;
      Memory::add<Data>(capacity);
   }

   /*
   * Deallocate a C array that is owned by this container.
   */
   template <typename Data>
   void FftwDRArray<Data>::deallocate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(!ref_.isAssociated());
      fftw_free(data_);
      data_ = nullptr;
      Memory::sub<Data>(capacity_);
      capacity_ = 0;
      UTIL_CHECK(!refCounter_.hasRefs());
   }

   /*
   * Associate this object with a slice of a different FftwDRArray.
   */
   template <typename Data>
   void FftwDRArray<Data>::associate(FftwDRArray<Data>& owner,
                                 int beginId, int capacity)
   {
      UTIL_CHECK(owner.isAllocated());
      UTIL_CHECK(owner.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= owner.capacity());
      UTIL_CHECK(!data_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      data_ = owner.cArray() + beginId;
      capacity_ = capacity;

      // Associate private ReferencecCounter of the data owner with the
      // CountedReference ref_ member variable of this data user.
      ref_.associate(owner.refCounter_);

      // On exit from CountedReference::associate, the ReferenceCounter
      // of the data owner is incremented and the CountedReference of
      // this data user has a pointer to that ReferenceCounter.
   }

   /*
   * Dissociate this object from array slice owned by another object.
   */
   template <typename Data>
   void FftwDRArray<Data>::dissociate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(ref_.isAssociated());

      data_ = nullptr;
      capacity_ = 0;
      ref_.dissociate(); // decrements counter mainained by owner
   }

} // namespace Cpu
} // namespace Prdc
} // namespace Pscf
#endif
