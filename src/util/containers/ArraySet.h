#ifndef UTIL_ARRAY_SET_H
#define UTIL_ARRAY_SET_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/containers/PArray.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   template <typename Data> class PArrayIterator;

   /**
   * A container for pointers to a subset of elements of an associated array.
   *
   * An ArraySet is a PArray that stores pointers to a subset of the
   * elements of an associated Array container or bare C array. Pointers 
   * to the elements of this set are stored in a contiguous sequence, with 
   * indices in the range 0, ..., size() - 1. The order in which these
   * pointers are stored is mutable, and generally changes whenever an 
   * element is removed.
   *
   * The append() method appends a pointer to a new element to the 
   * end of the sequence and increments the size. The remove() method 
   * removes a specified element, then moves the pointer of the last 
   * element to the space vacated by the removed element (unless the
   * removed element was the last in the sequence), and decrements the 
   * size. The order in which the remaining elements of an ArraySet 
   * are stored thus can change whenever an element is removed.
   *
   * An ArraySet provides O(N) sequential access to all elements of a
   * set, O(1) insertion and deletion, and O(1) access to a randomly 
   * chosen element.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class ArraySet : public PArray<Data>
   { 

      using PArray<Data>::ptrs_;
      using PArray<Data>::capacity_;
      using PArray<Data>::size_;

   public: 
   
      /**
      * Constructor.
      */
      ArraySet();
   
      /**
      * Destructor.
      */
      virtual ~ArraySet();
    
      /**
      * Associate with a C array and allocate required memory.
      *
      * This method associates an ArraySet with a bare C array, and 
      * allocates all memory required by the ArraySet. 
      *
      * An ArraySet may only be allocated once. This method throws an
      * Exception if it is called more than once.
      *
      * \param array    associated C array of Data objects
      * \param capacity number of elements in the array
      */
      void allocate(const Data* array, int capacity);
  
      /**
      * Associate with an Array container and allocate required memory.
      *
      * Invokes allocate(&array[0], array.capacity()) internally.
      *
      * \param array associated Array<Data> container
      */
      void allocate(const Array<Data>& array);
  
      /// \name Mutators
      //@{

      /**
      * Append an element to the set.
      *
      * This appends a new element to the end of the sequence.
      * This does not change the order of other elements.
      * 
      * \param data array element to be added.
      */
      void append(Data& data);

      /**
      * Remove an element from the set.
      *
      * Removal of an element generally changes the order of
      * the remaining elements.
      *
      * Throws an Exception if data is not in this ArraySet.
      *
      * \param data array element to be added.
      */
      void remove(const Data& data);

      /**
      * Pop the topmost from the set.
      *
      * Popping the top element does not change the order of
      * the remaining elements. 
      */
      Data& pop();

      /**
      * Reset to empty state.
      */ 
      void clear();

      //@}

      /// \name Accessors
      //@{
 
      /**
      * Return the current index of an element within the set, if any.
      *
      * This method returns the current index of the pointer to object
      * data within this ArraySet, in the range 0 < index < size() -1.
      * The method returns -1 if data is an element of the associated 
      * array but is not in the ArraySet.
      *
      * Throws an exception if data is not in the associated array.
      *
      * \param  data array element of interest.
      * \return current index of pointer to element within this ArraySet.
      */
      int index(const Data& data) const;

      /**
      * Return true if the ArraySet is initialized, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Return true if the ArraySet is valid, or throw an exception.
      */
      bool isValid() const;

      /**
      * Write the internal state of the ArraySet to std::cout
      */
      void dump() const;

      //@}

   private:

      // Associated C array of Data
      const Data* data_;
  
      // Locations of specific array elements within ptrs_
      int* tags_;
   
      // Return the array index in data_ of a Data* pointer
      int id(const Data *ptr) const;

      /// Copy constructor, declared private to prohibit copying.
      ArraySet(const ArraySet&);

      /// Assignment, declared private to prohibit assignment.
      ArraySet& operator = (const ArraySet&);

   }; 


   /* 
   * Default constructor.
   */
   template <typename Data>
   ArraySet<Data>::ArraySet()
    : data_(0),    
      tags_(0)
   {}

   /* 
   * Destructor.
   */
   template <typename Data>
   ArraySet<Data>::~ArraySet()
   {
      if (ptrs_) {
         Memory::deallocate<Data*>(ptrs_, capacity_);
      }
      if (tags_) {
         Memory::deallocate<int>(tags_, capacity_);
      }
   }
 
   /* 
   * Create an association with a C array, and allocate required memory.
   */
   template <typename Data>
   void ArraySet<Data>::allocate(const Data* array, int capacity) 
   {

      // Preconditions
      if (capacity == 0) UTIL_THROW("Zero capacity");
      if (capacity < 0)  UTIL_THROW("Negative capacity");
      if (array == 0) UTIL_THROW("Null array pointer");
      if (ptrs_) UTIL_THROW("ptrs_ array already allocated");
      if (tags_) UTIL_THROW("tags_ array already allocated");

      data_  = array;
      Memory::allocate<Data*>(ptrs_, capacity);
      Memory::allocate<int>(tags_, capacity);
      capacity_ = capacity;
      size_ = 0;

      for (int i = 0; i < capacity_; ++i) {
         ptrs_[i] =  0;
         tags_[i] = -1;
      }

   }

   /* 
   * Create association with an Array, and allocate required memory.
   */
   template <typename Data>
   void ArraySet<Data>::allocate(const Array<Data>& array) 
   {  allocate(&array[0], array.capacity()); }

   /*
   * Append an element to the set.
   */
   template <typename Data>
   void ArraySet<Data>::append(Data& data)
   {
      Data* const ptr = &data;

      // Check preconditions
      if (ptr < data_ || ptr >= data_ + capacity_) {
         UTIL_THROW("Pointer out of range");
      }
      int arrayId = id(ptr);
      int setId   = tags_[arrayId];
      if (setId >= 0) {
         UTIL_THROW("Attempt to add element that is already in set");
      }

      // Add to top of set
      ptrs_[size_]   = ptr;
      tags_[id(ptr)] = size_;
      ++size_;
   }

   /*
   * Remove a specific element from the set.
   */
   template <typename Data>
   void ArraySet<Data>::remove(const Data& data)
   {
      const Data* const ptr = &data;

      // Preconditions
      if (ptr < data_ || ptr >= data_ + capacity_) {
         UTIL_THROW("Pointer out of range");
      }
      int arrayId = id(ptr);
      int setId = tags_[arrayId];
      if (setId < 0) {
         UTIL_THROW("Element is not in set");
      }

      // Debug mode preconditions
      assert(setId < capacity_);
      assert(ptrs_[setId] == ptr);

      // Remove from set 
      tags_[arrayId] = -1;
      if (setId != size_ - 1) {
         ptrs_[setId] = ptrs_[size_-1];
         tags_[id(ptrs_[setId])] = setId;
      }
      ptrs_[size_-1] = 0;
      --size_;

   }

   /*
   * Remove the topmost element from the set.
   */
   template <typename Data>
   Data& ArraySet<Data>::pop()
   {
      // Precondition
      if (size_ == 0) {
         UTIL_THROW("Attempt to pop from empty ArraySet");
      }

      Data* ptr = ptrs_[size_ - 1];
      tags_[id(ptr)] = -1;
      ptrs_[size_-1] =  0;
      --size_;
      return *ptr;
   }

   /*
   * Reset to empty state.
   */
   template <typename Data>
   void ArraySet<Data>::clear() 
   { 
      assert(ptrs_ != 0);
      assert(tags_ != 0);
      size_ = 0;
      for (int i = 0; i < capacity_; ++i) {
         ptrs_[i] =  0;
         tags_[i] = -1;
      }
   }

   /**
   * Return the current index of an element within the set, or return
   * a negative value -1 if the element is not in the set.
   */
   template <typename Data>
   int ArraySet<Data>::index(const Data& data) const
   {
      const Data* const ptr = &data;
      if (ptr < data_ || ptr >= data_ + capacity_) {
         UTIL_THROW("Pointer out of range");
      }
      return tags_[id(ptr)];
   }

   /* 
   * Is this allocated?
   */
   template <typename Data>
   inline bool ArraySet<Data>::isAllocated() const
   { return ptrs_ != 0; }

   /* 
   * Return true if valid, or throw an exception.
   */
   template <typename Data>
   bool ArraySet<Data>::isValid() const
   {
      int i, j, size;

      if (ptrs_ != 0) {
         if (tags_ == 0) {
            UTIL_THROW("ptrs_ is allocated but tags_ is not");
         }

         // Loop over tags_
         size = 0;
         for (i = 0; i < capacity_ ; ++i) {
            j = tags_[i];
            if (j >= 0) {
               if (j >= capacity_) {
                  UTIL_THROW("tags_[i] > capacity_");
               }
               if (id(ptrs_[j]) != i) {
                  UTIL_THROW("Inconsistent tags and ptrs: 1");
               }
               ++size;
            }
         }
         if (size != size_) {
            UTIL_THROW("Number of nonnegative tags != size");
         }

         // Loop over ptrs in set
         size = 0;
         for (i = 0; i < size_ ; ++i) {
            if (ptrs_[i] !=0) {
               j = id(ptrs_[i]);
               if (tags_[j] != i) {
                  UTIL_THROW("Inconsistent tags and ptrs: 2");
               }
               ++size;
            } else {
               UTIL_THROW("Null pointer ptrs_[i] for i < size_");
            }
         }
         if (size != size_) {
            UTIL_THROW("Number of non-Null pointers != size");
         }

         // Check empty elements of ptrs_ array
         for (i = size_; i < capacity_ ; ++i) {
            if (ptrs_[i] != 0) {
               UTIL_THROW("Non-null pointer ptrs_[i] for i >= size_");
            }
         }

      } else {
         if (tags_ != 0) {
            UTIL_THROW("ptrs_ == 0, but tags_ != 0");
         }
         if (capacity_ != 0) {
            UTIL_THROW("ptrs_ == 0, but capacity_ != 0");
         }
         if (size_ != 0) {
            UTIL_THROW("ptrs_ == 0, but size_ != 0");
         }
      }
      return true;
   }

   // Private inline function

   template <typename Data>
   inline int ArraySet<Data>::id(const Data *ptr) const
   {  return int(ptr - data_); }

} 
#endif
