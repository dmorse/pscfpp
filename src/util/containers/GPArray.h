#ifndef UTIL_G_P_ARRAY_H
#define UTIL_G_P_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArray.h>
#include <util/global.h>

namespace Util
{

   /**
   * An automatically growable PArray.
   *
   * An GPArray is a PArray that grows as needed as objects are appended.
   * Like any PArray, it holds pointers to objects, rather than objects.
   * The associated objects are not destroyed when a PArray is deallocated
   * or destroyed.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class GPArray : public PArray<Data>
   {

   public:

      /**
      * Constructor.
      */
      GPArray();

      /**
      * Copy constructor, copy pointers.
      *
      * Allocates new Data* array and copies pointers to Data objects.
      *
      *\param other the GPArray to be copied.
      */
      GPArray(const GPArray<Data>& other);
   
      /**
      * Assignment, element by element.
      *
      * Preconditions: 
      * - Both this and other GPArrays must be allocated.
      * - Capacity of this GPArray must be >= size of RHS GPArray.
      *
      * \param other the rhs GPArray 
      */
      GPArray<Data>& operator=(const GPArray<Data>& other);

      /**
      * Destructor.
      *
      * Deletes array of pointers, if allocated previously.
      * Does not delete the associated Data objects.
      */
      virtual ~GPArray();

      /**
      * Append an element to the end of the sequence.
      *
      * Resizes array if space is inadequate. 
      *
      * \param data Data object to be appended
      */
      void append(Data& data);

      /**
      * Reserve memory for specified number of elements.
      *
      * Resizes and copies array if requested capacity is less than the
      * current capacity. Does nothing if requested capacity is greater
      * than current capacity.
      *
      * \param capacity number of elements for which to reserve space.
      */
      void reserve(int capacity);

      /**
      * Deallocate (delete) underlying array of pointers.
      */
      void deallocate();

      /**
      * Reset to empty state.
      */ 
      void clear();

      /**
      * Is this GPArray allocated?
      */
      bool isAllocated() const;

   protected:

      using PArray<Data>::ptrs_;
      using PArray<Data>::capacity_;
      using PArray<Data>::size_;

   };

   /*
   * Default constructor.
   */
   template <typename Data>
   inline GPArray<Data>::GPArray()
    : PArray<Data>()
   {}

   /**
   * Copy constructor, copy pointers.
   *
   * Allocates a new Data* array and copies all pointer values.
   *
   *\param other the GPArray to be copied.
   */
   template <typename Data>
   GPArray<Data>::GPArray(const GPArray<Data>& other) 
    : PArray<Data>()
   {
      assert(other.capacity_ >= other.size_);
      if (other.ptrs_ == 0) {
         assert(other.capacity_ == 0);
         assert(other.size_ == 0);
         ptrs_ = 0;
         capacity_ = 0;
         size_ = 0;
      } else { 
         assert(other.capacity_ > 0);
         // Allocate array of Data* pointers
         Memory::allocate<Data*>(ptrs_, other.capacity_);
         capacity_ = other.capacity_;
         size_ = other.size_;
         // Copy pointers
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               ptrs_[i] = other.ptrs_[i];
            }
         }
         // Nullify unused elements of ptrs_ array
         if (capacity_ > size_) {
            for (int i = size_; i < capacity_; ++i) {
               ptrs_[i] = 0;
            }
         }
      }
   }

   /*
   * Assignment, element by element.
   */
   template <typename Data>
   GPArray<Data>& GPArray<Data>::operator=(const GPArray<Data>& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      clear();
      for (int i = 0; i < other.size_; ++i) {
         append(other[i]);
      }
      return *this;
   }

   /*
   * Destructor.
   */
   template <typename Data>
   GPArray<Data>::~GPArray()
   {
      size_ = 0; 
      if (isAllocated()) {
         Memory::deallocate(ptrs_, capacity_);
         capacity_ = 0; 
      }
   }

   /*
   * Reserve space for the underlying array of Data* pointers.
   */
   template <typename Data>
   void GPArray<Data>::reserve(int capacity) 
   {
      if (capacity <= 0) {
         UTIL_THROW("Cannot reserve with capacity <=0");
      }
      if (ptrs_ == 0) {
         assert(capacity_ == 0);
         assert(size_ == 0);
         Memory::allocate<Data*>(ptrs_, capacity);
         capacity_ = capacity;
         size_ = 0;
      } else if (capacity > capacity_) {
         assert(capacity_ > 0);
         assert(capacity_ >= size_);
         assert(size_ >= 0);
         Data** newPtr = 0;
         Memory::allocate<Data*>(newPtr, capacity);
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               newPtr[i] = ptrs_[i];
            }
         }
         Memory::deallocate<Data*>(ptrs_, capacity_);
         ptrs_ = newPtr;
         capacity_ = capacity;
      }
   }

   /*
   * Deallocate associated memory.
   */
   template <typename Data>
   void GPArray<Data>::deallocate() 
   {  
      size_ = 0; 
      if (isAllocated()) {
         assert(capacity_ > 0);
         Memory::deallocate<Data*>(ptrs_, capacity_);
         capacity_ = 0; 
      } 
   }

   /*
   * Append an element to the end of the PArray.
   */
   template <typename Data>
   void GPArray<Data>::append(Data& data) 
   {
      if (!isAllocated()) {
         assert(capacity_ == 0);
         Memory::allocate<Data*>(ptrs_, 64);
         capacity_ = 64;
         size_ = 0;
      } else if (size_ == capacity_) {
         assert(capacity_ > 0);
         assert(capacity_ >= size_);
         assert(size_ >= 0);
         Data** newPtr = 0;
         Memory::allocate<Data*>(newPtr, 2*capacity_);
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               newPtr[i] = ptrs_[i];
            }
         }
         Memory::deallocate<Data*>(ptrs_, capacity_);
         ptrs_ = newPtr;
         capacity_ = 2*capacity_;
         // size_ is unchanged
      }
      // Append new element
      ptrs_[size_] = &data;
      ++size_;
      assert(size_ <= capacity_);
   }

   /*
   * Reset to empty state, without deallocating.
   */
   template <typename Data>
   inline void GPArray<Data>::clear()
   {  size_ = 0; }

   /*
   * Is this GPArray allocated?
   */
   template <class Data>
   inline bool GPArray<Data>::isAllocated() const
   {  return (bool)ptrs_; }

} 
#endif
