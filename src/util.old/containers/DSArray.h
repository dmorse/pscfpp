#ifndef UTIL_DS_ARRAY_H
#define UTIL_DS_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ArrayIterator.h>
#include <util/containers/ConstArrayIterator.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * Dynamically allocated array with variable logical size.
   *
   * A DSArray < Data > is a wrapper for a dynamically allocated C array,
   * with continuous elements and a logical size that may be less than
   * or equal to its physical capacity. The logical size is the number
   * of contiguous elements that have been added using the append() method.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class DSArray
   {

   public:

      /**
      * Constructor.
      */
      DSArray();

      /**
      * Copy constructor.
      *
      *\param other the DSArray to be copied.
      */
      DSArray(const DSArray< Data >& other);

      /**
      * Assignment, element by element.
      *
      * Capacity of LHS must be either zero or equal that of RHS DSArray.
      *
      * \param other the RHS DSArray
      */
      DSArray<Data>& operator=(const DSArray<Data>& other);

      /**
      * Destructor.
      */
      virtual ~DSArray();

      /**
      * Allocates the underlying C array.
      *
      * Throw an exception if the DSArray has already been
      * allocated - A DSArray can only be allocated once.
      *
      * \param capacity number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Append data to the end of the array.
      *
      * \param data Data to add to end of array.
      */
      void append(const Data &data);

      /**
      * Modify logical size without modifying data.
      *
      * The size parameter must be non-negative and may not exceed
      * the physical allocated capacity.
      *
      * This function simply changes the logical size without
      * modifying any elements of the underlying physical array.
      * When the size increases, added elements are uninitialized.
      *
      * \param size new logical size, 0 <= size < capacity.
      */
      void resize(int size);

      /**
      * Set logical size to zero.
      */
      void clear();

      /**
      * Serialize a DSArray to/from an Archive.
      *
      * \param ar        archive
      * \param version   archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Set an ArrayIterator to the beginning of this Array.
      *
      * \param iterator ArrayIterator, initialized on output.
      */
      void begin(ArrayIterator<Data> &iterator);

      /**
      * Set a ConstArrayIterator to the beginning of this Array.
      *
      * \param iterator ConstArrayIterator, initialized on output.
      */
      void begin(ConstArrayIterator<Data> &iterator) const;

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

      /**
      * Return physical capacity of array.
      */
      int capacity() const;

      /**
      * Return logical size of this array (i.e., number of elements).
      */
      int size() const;

      /**
      * Return true if the DSArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

   protected:

      /// C array of Data elements.
      Data *data_;

      /// Logical size of array (number of elements used).
      int size_;

      /// Capacity (physical size) of underlying C array.
      int capacity_;

   };

   // Method definitions

   /*
   * Constructor.
   */
   template <class Data>
   DSArray<Data>::DSArray()
    : data_(0),
      size_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor.
   */
   template <class Data>
   DSArray<Data>::DSArray(const DSArray< Data >& other)
    : data_(0),
      size_(0),
      capacity_(0)
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DSArray must be allocated.");
       }

      Memory::allocate<Data>(data_, other.capacity_);
      capacity_ = other.capacity_;
      size_ = other.size_;
      for (int i = 0; i < size_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment, element by element.
   *
   * Capacity of LHS DSArray must be zero or == capacity of RHS DSArray.
   */
   template <class Data>
   DSArray<Data>& DSArray<Data>::operator=(const DSArray<Data>& other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DSArray must be allocated.");
      }

      if (!isAllocated()) {
         allocate(other.capacity_);
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign DSArrays of unequal capacity");
      }

      // Copy elements and set size
      for (int i = 0; i < other.size_; ++i) {
         data_[i] = other[i];
      }
      size_ = other.size_;

      return *this;
   }

   /*
   * Destructor.
   */
   template <class Data>
   DSArray<Data>::~DSArray()
   {
       size_ = 0;
       if (isAllocated()) {
          assert(capacity_);
          Memory::deallocate<Data>(data_, capacity_);
          capacity_ = 0;
       }
   }

   /*
   * Allocates the underlying C array.
   */
   template <class Data>
   void DSArray<Data>::allocate(int capacity)
   {
      if (isAllocated()) {
         UTIL_THROW("Cannot re-allocate a DSArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate a DSArray with capacity <= 0");
      }
      Memory::allocate<Data>(data_, capacity);
      capacity_ = capacity;
      size_ = 0;
   }

   /*
   * Serialize a DSArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void DSArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         if (!isAllocated()) {
            UTIL_THROW("Cannot save unallocated DSArray");
         }
         capacity = capacity_;
      }
      ar & capacity;
      ar & size_;
      if (Archive::is_loading()) {
         if (capacity <= 0) {
            UTIL_THROW("Invalid DSArray input capacity on load, capacity <= 0");
         }
         if (size_ < 0) {
            UTIL_THROW("Invalid DSArray input size on load, size < 0");
         }
         if (!isAllocated()) {
            allocate(capacity);
         }
         if (size_ > capacity_) {
            UTIL_THROW("Inconsistent DSArray size and capacity on load");
         }
      }
      for (int i = 0; i < size_; ++i) {
         ar & data_[i];
      }
   }

   /**
   * Set an ArrayIterator to the beginning of this Array.
   *
   * \param iterator ArrayIterator, initialized on output.
   */
   template <class Data>
   inline
   void DSArray<Data>::begin(ArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Set a ConstArrayIterator to the beginning of this Array.
   */
   template <class Data>
   inline
   void DSArray<Data>::begin(ConstArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline Data& DSArray<Data>::operator[] (int i)
   {
      assert(i < size_);
      assert(i >= 0);
      return data_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline const Data& DSArray<Data>::operator[] (int i) const
   {
      assert(i < size_);
      assert(i >= 0 );
      return data_[i];
   }

   /*
   * Append data to the end of the array.
   */
   template <class Data>
   inline void DSArray<Data>::append(const Data &data)
   {
      if (size_ == capacity_) {
         UTIL_THROW("Attempt to add to full DSArray");
      }
      data_[size_] = data;
      ++size_;
   }

   /**
   * Modify logical size without modifying data.
   *
   * The size parameter must be non-negative and may not exceed
   * the capacity.
   *
   * This function simply changes the logical size of without
   * modifying any elements of the underlying physical array.
   * If the size increases, added elements are uninitialized.
   *
   * \param size new logical size, 0 <= size < capacity.
   */
   template <class Data>
   inline void DSArray<Data>::resize(int size)
   {  size_ = size; }

   /*
   * Set logical size to zero.
   */
   template <class Data>
   inline void DSArray<Data>::clear()
   {  size_ = 0; }

   /*
   * Return physical capacity of array.
   */
   template <class Data>
   inline int DSArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return logical size of this array (i.e., number of elements).
   */
   template <class Data>
   inline int DSArray<Data>::size() const
   {  return size_; }

   /*
   * Return true if the DSArray has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DSArray<Data>::isAllocated() const
   {  return (bool)data_; }

}
#endif
