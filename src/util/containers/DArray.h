#ifndef UTIL_D_ARRAY_H
#define UTIL_D_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * Dynamically allocatable contiguous array template.
   *
   * A DArray wraps a dynamically allocated C Array. Unlike an STL
   * std::vector, a DArray cannot be resized after it is allocated.
   * Addresses of elements thus remain constant until de-allocation.
   *
   * The Array<Data> base class provides bounds checking when compiled
   * in debug mode. 
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class DArray : public Array<Data>
   {

      using Array<Data>::data_;
      using Array<Data>::capacity_;

   public:

      /**
      * Default constructor.
      */
      DArray();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the DArray to be copied.
      */
      DArray(const DArray<Data>& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DArray();

      /**
      * Assignment operator.
      *
      * If this DArray is not allocated, allocates and copies all elements.
      *
      * If this and the other DArray are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS DArray
      */
      DArray<Data>& operator = (const DArray<Data>& other);

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the DArray is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the DArray is not allocated.
      */
      void deallocate();

      /**
      * Return true if the DArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Serialize a DArray to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   };


   /*
   * Constructor.
   */
   template <class Data>
   DArray<Data>::DArray()
    : Array<Data>()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the DArray to be copied.
   */
   template <class Data>
   DArray<Data>::DArray(const DArray<Data>& other)
    : Array<Data>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other DArray must be allocated.");
      }
      Memory::allocate(data_, other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Destructor.
   */
   template <class Data>
   DArray<Data>::~DArray()
   {
      if (isAllocated()) {
         Memory::deallocate<Data>(data_, capacity_);
         capacity_ = 0;
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other DArray is not allocated.
   * \throw Exception if both DArrays are allocated with unequal capacities.
   *
   * \param other the rhs DArray
   */
   template <class Data>
   DArray<Data>& DArray<Data>::operator = (const DArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DArray must be allocated.");
      }

      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign DArrays of unequal capacity");
      }

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the DArray has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <class Data>
   void DArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a DArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      Memory::allocate<Data>(data_, capacity);
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this DArray is not allocated.
   */
   template <class Data>
   void DArray<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      Memory::deallocate<Data>(data_, capacity_);
      capacity_ = 0;
   }

   /*
   * Return true if the DArray has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DArray<Data>::isAllocated() const
   {  return (bool)data_; }

   /*
   * Serialize a DArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void DArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent DArray capacities");
            }
         }
      }
      if (isAllocated()) {
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

}
#endif
