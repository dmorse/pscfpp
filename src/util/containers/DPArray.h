#ifndef UTIL_D_P_ARRAY_H
#define UTIL_D_P_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArray.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * A dynamic array that only holds pointers to its elements.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class DPArray : public PArray<Data>
   {

   public:

      /**
      * Constructor.
      */
      DPArray();

      /**
      * Copy constructor, copy pointers.
      *
      * Allocates new Data* array and copies pointers to Data objects.
      *
      *\param other the DPArray to be copied.
      */
      DPArray(const DPArray<Data>& other);
   
      /**
      * Destructor.
      *
      * Deletes array of pointers, if allocated previously.
      * Does not delete the associated Data objects.
      */
      virtual ~DPArray();

      /**
      * Assignment, element by element.
      *
      * Preconditions: 
      * - Both this and other DPArrays must be allocated.
      * - Capacity of this DPArray must be >= size of RHS DPArray.
      *
      * \param other the rhs DPArray 
      */
      DPArray<Data>& operator=(const DPArray<Data>& other);

      /**
      * Allocate an array of pointers to Data.
      *
      * Throw an Exception if the DPArray has already been
      * allocated - A DPArray can only be allocated once.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Append an element to the end of the sequence.
      *
      * \param data Data object to be appended
      */
      void append(Data& data);

      /**
      * Reset to empty state.
      */ 
      void clear();

      /**
      * Is this DPArray allocated?
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
   inline DPArray<Data>::DPArray()
    : PArray<Data>()
   {}

   /**
   * Copy constructor, copy pointers.
   *
   * Allocates a new Data* array and copies all pointer values.
   *
   *\param other the DPArray to be copied.
   */
   template <typename Data>
   DPArray<Data>::DPArray(const DPArray<Data>& other) 
    : PArray<Data>()
   {
      if (other.size_ > other.capacity_) {
         UTIL_THROW("Inconsistent size and capacity");
      }
      if (other.isAllocated()) {
         // Allocate array of Data* pointers
         Memory::allocate<Data*>(ptrs_, other.capacity_);
         capacity_ = other.capacity_;
         // Copy pointers
         for (int i = 0; i < other.size_; ++i) {
            ptrs_[i] = other.ptrs_[i];
         }
         size_ = other.size_;
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
   DPArray<Data>& DPArray<Data>::operator=(const DPArray<Data>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      // Preconditions
      if (ptrs_ == 0) {
         UTIL_THROW("LHS DPArray in assignment is not allocated");
      }
      if (other.ptrs_ == 0) {
         UTIL_THROW("RHS DPArray in assignment is not allocated");
      }
      if (capacity_ < other.size_) {
         UTIL_THROW("LHS DPArray is too small");
      }

      // Copy pointers
      int i;
      for (i = 0; i < other.size_; ++i) {
         ptrs_[i] = other[i];
      }
      size_ = other.size_;

      // Nullify any unused elements
      if (capacity_ > size_) {
         for (i = size_; i < capacity_; ++i) {
            ptrs_[i] = 0;
         }
      }

      return *this;
   }

   /*
   * Destructor.
   */
   template <typename Data>
   DPArray<Data>::~DPArray()
   {
      size_ = 0;
      if (ptrs_) {
         assert(capacity_ > 0);
         Memory::deallocate<Data*>(ptrs_, capacity_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying array of Data* pointers.
   */
   template <typename Data>
   void DPArray<Data>::allocate(int capacity) 
   {
      // Preconditions
      if (!(ptrs_ == 0)) {
         UTIL_THROW("Cannot re-allocate a DPArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate a DPArray with capacity <=0");
      }
      Memory::allocate<Data*>(ptrs_, capacity);
      capacity_ = capacity;
   }

   /*
   * Append an element to the end of the PArray.
   */
   template <typename Data>
   inline void DPArray<Data>::append(Data& data) 
   {
      if (!isAllocated()) {
         UTIL_THROW("Error: Attempt to append to unallocated DPArray");
      }
      if (size_ == capacity_) {
         UTIL_THROW("Error: Attempt to append data to a full DPArray");
      }
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Reset to empty state.
   */
   template <typename Data>
   void DPArray<Data>::clear()
   {
      assert(isAllocated());
      size_ = 0;
   }

   /*
   * Is this DPArray allocated?
   */ 
   template <typename Data>
   inline 
   bool DPArray<Data>::isAllocated() const
   { return (bool)ptrs_;  }



} 
#endif
