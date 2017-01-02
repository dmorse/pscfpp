#ifndef UTIL_S_SET_H
#define UTIL_S_SET_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArrayIterator.h>
#include <util/containers/ConstPArrayIterator.h>
#include <util/global.h>

namespace Util
{

   /**
   * Statically allocated array of pointers to an unordered set.
   *
   * An SSet is a statically allocated array that holds pointers to 
   * a set of objects. It implements the same interface as PArray 
   * and FPArray, plus additional remove() and index() methods. As
   * for any pointer array container, the [] operator returns an 
   * associated object by reference .
   *
   * An SSet holds a set of pointers in a contiguous array. The size
   * is the number of pointers now in the container, and the Capacity
   * is the maximum number it can hold. The class is implemented as a
   * wrapper for a statically allocated C array of Capacity elements.
   *
   * The append method adds a pointer to the end of the sequence. 
   * The remove method removes an object from the set, or throws an 
   * exception if the object is not found in the set. As for an 
   * ArraySet, the remove method repacks the sequence of pointers by 
   * moving the last element to the position of the element that is
   * being removed. Removal of an element thus generally changes the 
   * order in which the remaining elements are stored.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data, int Capacity>
   class SSet 
   {

   public:

      /**
      * Default constructor.
      */
      SSet();

      /**
      * Copy constructor.
      *
      * Copies all pointers.
      *
      *\param other the SSet to be copied.
      */
      SSet(const SSet<Data, Capacity>& other);

      /**
      * Assignment, element by element.
      *
      * \param other the rhs SSet 
      */
      SSet<Data, Capacity>& operator=(const SSet<Data, Capacity>& other);

      /**
      * Destructor.
      */
      ~SSet();

      /**
      * Add an object to the set.
      *
      * Appends a pointer to the object to the end of the sequence.
      *
      * \param data Data to add to end of array.
      */
      void append(Data &data);

      /**
      * Remove an object from the set.
      *
      * Removal of an object generally changes the storage order
      * of the remaining objects.
      *
      * \throw Exception if object data is not in the Set.
      *
      * \param data object to be removed.
      */
      void remove(const Data& data);

      /**
      * Set logical size to zero and nullify all elements.
      */
      void clear();

      /**
      * Return physical capacity of array.
      */
      int capacity() const;

      /**
      * Return logical size of this array.
      */
      int size() const;

      /**
      * Is an object an element of the set?
      *
      * \param  data object of interest.
      */
      bool isElement(const Data& data) const;

      /**
      * Return the current index of an object within the set, if any.
      *
      * This method returns the current index of the pointer to object
      * data within this SSet, in the range 0 < index < size() - 1.
      * The method returns -1 if data is the object is not in the set.
      *
      * Throws an exception if data is not in the associated array.
      *
      * \param  data object of interest.
      * \return current index of pointer to element within this SSet.
      */
      int index(const Data& data) const;

      /**
      * Set a PArrayIterator to the beginning of this Array.
      *
      * \param iterator PArrayIterator, initialized on output. 
      */
      void begin(PArrayIterator<Data> &iterator);

      /**
      * Set a ConstPArrayIterator to the beginning of this Array.
      *
      * \param iterator ConstPArrayIterator, initialized on output. 
      */
      void begin(ConstPArrayIterator<Data> &iterator) const;

      /**
      * Mimic C array subscripting.
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i);

      /**
      * Mimic C array subscripting.
      *
      * \param i array index
      * \return const reference to element i
      */
      const Data& operator[] (int i) const;

   protected:

      /// Array of pointers to Data objects.
      Data* ptrs_[Capacity];

      /// Logical size of array (number of elements in array).
      int  size_;

   };
      
   // Method definitions

   /*
   * Constructor.
   */
   template <typename Data, int Capacity>
   inline SSet<Data, Capacity>::SSet()
    : size_(0)
   {}

   /*
   * Copy constructor, copy all pointers.
   */
   template<typename Data, int Capacity>
   SSet<Data, Capacity>::SSet(const SSet<Data, Capacity>& other) 
   {

      // Copy pointer values
      int i;
      for (i = 0; i < other.size_; ++i) {
         ptrs_[i] = other.ptrs_[i];
      }
      size_ = other.size_;

      // Nullify any unused elements
      if (Capacity > size_) {
         for (i = size_; i < Capacity; ++i) {
            ptrs_[i] = 0;
         }
      }

   }

   /*
   * Assignment, element by element.
   */
   template <typename Data, int Capacity>
   SSet<Data, Capacity>& 
   SSet<Data, Capacity>::operator=(const SSet<Data, Capacity>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (Capacity < other.size_) {
         UTIL_THROW("LHS SSet is too small");
      }

      // Copy pointers
      for (int i = 0; i < other.size_; ++i) {
         ptrs_[i] = other[i];
      }
      size_ = other.size_;

      // Nullify any unused elements
      if (Capacity > size_) {
         for (int i = size_; i < Capacity; ++i) {
            ptrs_[i] = 0;
         }
      }

      return *this;
   }
 
   /*
   * Destructor.
   */
   template <typename Data, int Capacity>
   SSet<Data, Capacity>::~SSet()
   {}

   /*
   * Return physical capacity of array.
   */
   template <typename Data, int Capacity>
   inline int SSet<Data, Capacity>::capacity() const
   { return Capacity; }

   /*
   * Return logical size of this array.
   */
   template <typename Data, int Capacity>
   inline int SSet<Data, Capacity>::size() const
   { return size_; }

   /*
   * Set a PArrayIterator to the beginning of this Array.
   *
   * \param iterator PArrayIterator, initialized on output. 
   */
   template <typename Data, int Capacity>
   inline void SSet<Data, Capacity>::begin(PArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(const_cast<Data**>(ptrs_));
      iterator.setEnd(const_cast<Data**>(ptrs_ + size_));
   }

   /*
   * Set a PArrayIterator to the beginning of this Array.
   */
   template <typename Data, int Capacity>
   inline void SSet<Data, Capacity>::begin(ConstPArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(const_cast<Data**>(ptrs_));
      iterator.setEnd(const_cast<Data**>(ptrs_ + size_));
   }

   /*
   * Mimic C array subscripting.
   */
   template <typename Data, int Capacity>
   inline Data& SSet<Data, Capacity>::operator[] (int i)
   {
      assert(i < size_);
      assert(i >= 0);
      return *ptrs_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <typename Data, int Capacity>
   inline const Data& SSet<Data, Capacity>::operator[] (int i) const
   {
      assert(i < size_);
      assert(i >= 0 );
      return *ptrs_[i];
   }

   /*
   * Append data to the end of the array.
   *
   * \param data Data to add to end of array.
   */
   template <typename Data, int Capacity>
   inline void SSet<Data, Capacity>::append(Data &data) 
   {
      if (size_ == Capacity) {
         UTIL_THROW("Attempt to add to full SSet");
      }
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Clear out the container.
   */
   template <typename Data, int Capacity>
   inline void SSet<Data, Capacity>::clear() 
   { 
      size_ = 0; 
      for (int i=0; i < Capacity; ++i) {
         ptrs_[i] = 0;
      }
   }
 
   /*
   * Remove an element from the set.
   */
   template <typename Data, int Capacity>
   void SSet<Data, Capacity>::remove(const Data& data)
   {
      const Data* const ptr = &data;
      int   i;
      bool  found;

      // Search for object
      found = false;
      for (i = 0; i < size_; ++i) {
         if (ptr == ptrs_[i]) {
            found = true;
            break;
         }
      }

      if (!found) {
         UTIL_THROW("Element is not in set");
      }

      // Remove element
      if (i != size_ - 1) {
         ptrs_[i] = ptrs_[size_-1];
      }
      ptrs_[size_ - 1] = 0;
      --size_;

   }

   /*
   * Return true if an object is in the set, false otherwise.
   */
   template <typename Data, int Capacity>
   bool SSet<Data, Capacity>::isElement(const Data& data) const
   {
      const Data* const ptr = &data;
      bool  found;

      // Search for element in set
      found = false;
      for (int i = 0; i < size_; ++i) {
         if (ptr == ptrs_[i]) {
            found = true;
            break;
         }
      }
      return found;
   }

   /**
   * Return the current index of an element within the set, 
   * or return -1 if the element is not in the set.
   */
   template <typename Data, int Capacity>
   int SSet<Data, Capacity>::index(const Data& data) const
   {
      const Data* const ptr = &data;
      int   i;
      bool  found;

      // Search for element in set
      found = false;
      for (i = 0; i < size_; ++i) {
         if (ptr == ptrs_[i]) {
            found = true;
            break;
         }
      }
      if (found) {
         return i;
      } else {
         return -1;
      }
   }

} 
#endif
