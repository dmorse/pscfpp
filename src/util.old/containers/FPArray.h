#ifndef UTIL_FP_ARRAY_H
#define UTIL_FP_ARRAY_H

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
   * Statically allocated pointer array.
   *
   * A FPArray is a statically allocated array that actually holds pointers
   * to Data objects, but for which the [] operator returns a reference to 
   * the associated object. It is implemented as a wrapper for a statically 
   * allocated C array of Data* pointers. A FPArray is not responsible for
   * destroying the associated Data objects.
   *
   * The interface of an FPArray is identical to that of an FSArray.
   * An FPArray has both a capacity that is set at compile time, which is 
   * the physical size of the underlying C array, and a size, which is the 
   * number of contiguous elements (indexed from 0 to size-1) that contain 
   * valid pointers.  The size can only be increased only by the append() 
   * method, which adds an element to the end of the array. 
   *
   * When compiled in debug mode, the operator [] checks that the index is 
   * less than the size and non-negative.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data, int Capacity>
   class FPArray 
   {

   public:

      /**
      * Default constructor.
      */
      FPArray();

      /**
      * Copy constructor.
      *
      * Copies all pointers.
      *
      *\param other the FPArray to be copied.
      */
      FPArray(const FPArray<Data, Capacity>& other);

      /**
      * Destructor.
      */
      ~FPArray();

      /**
      * Assignment, element by element.
      *
      * \param other the rhs FPArray 
      */
      FPArray<Data, Capacity>& operator=(const FPArray<Data, Capacity>& other);

      /**
      * Append an element to the end of the array.
      *
      * \param data Data to add to end of array.
      */
      void append(Data &data);

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
      * Set an iterator to begin this container.
      *
      * \param iterator PArrayIterator, initialized on output. 
      */
      void begin(PArrayIterator<Data> &iterator);

      /**
      * Set a const iterator to begin this container.
      *
      * \param iterator ConstPArrayIterator, initialized on output. 
      */
      void begin(ConstPArrayIterator<Data> &iterator) const;

      /**
      * Get an element by reference (mimic C-array subscripting).
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i);

      /**
      * Get an element by const reference (mimic C-array subscripting).
      *
      * \param i array index
      * \return const reference to element i
      */
      const Data& operator[] (int i) const;

   protected:

      /// Array of pointers to Data objects.
      Data* ptrs_[Capacity];

      /// Logical size of array (number of elements used).
      int  size_;

   };
      
   // Method definitions

   /**
   * Constructor.
   */
   template <typename Data, int Capacity>
   inline FPArray<Data, Capacity>::FPArray()
    : size_(0)
   {}

   /*
   * Copy constructor, copy all pointers.
   */
   template<typename Data, int Capacity>
   FPArray<Data, Capacity>::FPArray(const FPArray<Data, Capacity>& other) 
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
   FPArray<Data, Capacity>& 
   FPArray<Data, Capacity>::operator=(const FPArray<Data, Capacity>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (Capacity < other.size_) {
         UTIL_THROW("LHS FPArray is too small");
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
 
   /**
   * Destructor.
   */
   template <typename Data, int Capacity>
   FPArray<Data, Capacity>::~FPArray()
   {}

   /**
   * Return physical capacity of array.
   */
   template <typename Data, int Capacity>
   inline int FPArray<Data, Capacity>::capacity() const
   { return Capacity; }

   /*
   * Return logical size of this array.
   */
   template <typename Data, int Capacity>
   inline int FPArray<Data, Capacity>::size() const
   { return size_; }

   /*
   * Set a PArrayIterator to the beginning of this Array.
   */
   template <typename Data, int Capacity>
   inline 
   void FPArray<Data, Capacity>::begin(PArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(const_cast<Data**>(ptrs_));
      iterator.setEnd(const_cast<Data**>(ptrs_ + size_));
   }

   /*
   * Set a PArrayIterator to the beginning of this Array.
   */
   template <typename Data, int Capacity>
   inline 
   void FPArray<Data, Capacity>::begin(ConstPArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(const_cast<Data**>(ptrs_));
      iterator.setEnd(const_cast<Data**>(ptrs_ + size_));
   }

   /*
   * Mimic C array subscripting.
   *
   * \param  i array index
   * \return reference to element i
   */
   template <typename Data, int Capacity>
   inline Data& FPArray<Data, Capacity>::operator[] (int i)
   {
      assert(i < size_);
      assert(i >= 0);
      return *ptrs_[i];
   }

   /*
   * Mimic C array subscripting.
   *
   * \param i array index
   * \return const reference to element i
   */
   template <typename Data, int Capacity>
   inline const Data& FPArray<Data, Capacity>::operator[] (int i) const
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
   inline void FPArray<Data, Capacity>::append(Data &data) 
   {
      if (size_ == Capacity) {
         UTIL_THROW("Attempt to add to full FPArray");
      }
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Set logical size to zero.
   */
   template <typename Data, int Capacity>
   inline void FPArray<Data, Capacity>::clear() 
   { 
      size_ = 0; 
      for (int i=0; i < Capacity; ++i) {
         ptrs_[i] = 0;
      }
   }

} 
#endif
