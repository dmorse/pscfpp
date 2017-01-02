#ifndef UTIL_G_ARRAY_H
#define UTIL_G_ARRAY_H

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
   * An automatically growable array, analogous to a std::vector.
   *
   * An GArray is an array that grows as needed as objects are appended.
   * It has a logical size that grows when objects are appended, which 
   * is always less than or equal to the current physical capacity. If
   * an object is added when the size is already equal to the capacity,
   * the array will be resized and copied to a new location in memory.
   * The elements of a GArray are deleted when the GArray is destroyed
   * or deallocated.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class GArray 
   {

   public:

      /**
      * Constructor.
      */
      GArray();

      /**
      * Copy constructor, copy pointers.
      *
      * Allocates new C-array and copies pointers to Data objects.
      *
      *\param other the GArray to be copied.
      */
      GArray(const GArray<Data>& other);
   
      /**
      * Assignment, element by element.
      *
      * \param other the rhs GArray 
      */
      GArray<Data>& operator = (const GArray<Data>& other);

      /**
      * Destructor.
      *
      * Deletes array of pointers, if allocated previously.
      * Does not delete the associated Data objects.
      */
      virtual ~GArray();

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
      *
      * Sets capacity and size to zero.
      */
      void deallocate();

      /**
      * Reset to empty state.
      *
      * Sets size to zero, but leaves capacity unchanged.
      * Does not call destructor for deleted elements.
      */ 
      void clear();

      /**
      * Serialize a GArray to/from an Archive.
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
      * Append an element to the end of the sequence.
      *
      * Resizes array if space is inadequate. 
      *
      * \param data Data object to be appended
      */
      void append(const Data& data);

      /**
      * Resizes array so that it contains n elements.
      *
      * This function changes the size of the array to n, and changes 
      * the capacity iff necesary to accomodate the change in size. 
      * Upon return, size is set to n. In what follows, "size" and 
      * "capacity" refer to values on entry:
      *
      * If n < size, size is reset, but no destructors are called
      * If n > size, all added elements are value initialized
      * If n > capacity, new memory is allocated and the array is moved 
      *
      * \param n desired number of elements
      */
      void resize(int n);

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
      * Return logical size of this array (i.e., current number of elements).
      */
      int size() const;

      /**
      * Is this array allocated?
      */
      bool isAllocated() const;

   private:

      /// Array of Data elements.
      Data *data_;

      /// Logical size of array (number of elements used).
      int  size_;

      /// Maxium size of array
      int capacity_;

   }; // class GArray

   // Method definitions

   /*
   * Default constructor.
   */
   template <typename Data>
   GArray<Data>::GArray()
    : data_(0),
      size_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor, copies array elements.
   *
   * Allocates a new C-array and copies all elements.
   */
   template <typename Data>
   GArray<Data>::GArray(const GArray<Data>& other) 
    : data_(0),
      size_(0),
      capacity_(0)
   {
      assert(other.size_ <= other.capacity_);
      if (other.isAllocated()) {
         assert(other.capacity_ > 0);
         // Allocate new array
         Memory::allocate<Data>(data_, other.capacity_);
         capacity_ = other.capacity_;
         // Copy objects
         for (int i = 0; i < other.size_; ++i) {
            data_[i] = other.data_[i];
            ++size_;
         }
      }
   }

   /*
   * Destructor.
   */
   template <typename Data>
   GArray<Data>::~GArray()
   {
      size_ = 0;
      if (isAllocated()) {
         Memory::deallocate<Data>(data_, capacity_);
         capacity_ = 0;
      }
   }

   /*
   * Assignment, element by element.
   */
   template <typename Data>
   GArray<Data>& GArray<Data>::operator = (const GArray<Data>& other) 
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
   * Reserve space for the underlying array.
   */
   template <typename Data>
   void GArray<Data>::reserve(int capacity) 
   {
      if (capacity <= 0) {
         UTIL_THROW("Cannot reserve with capacity <=0");
      }
      if (!isAllocated()) {
         assert(capacity_ == 0);
         assert(size_ == 0);
         Memory::allocate<Data>(data_, capacity);
         capacity_ = capacity;
         size_ = 0;
      } else if (capacity > capacity_) {
         assert(capacity_ > 0);
         assert(capacity_ >= size_);
         Data* newPtr = 0;
         Memory::allocate<Data>(newPtr, capacity);
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               newPtr[i] = data_[i];
            }
         }
         Memory::deallocate<Data>(data_, capacity_);
         data_ = newPtr;
         capacity_ = capacity;
      }
   }

   /*
   * Delete associated C array.
   */
   template <typename Data>
   void GArray<Data>::deallocate() 
   {
      size_ = 0;
      if (isAllocated()) {
         Memory::deallocate<Data>(data_, capacity_);
         capacity_ = 0; 
      }
   }

   /*
   * Reset to empty state, without deallocating.
   */
   template <typename Data>
   void GArray<Data>::clear()
   {  size_ = 0; }

   /*
   * Append an element to the end of the Array.
   */
   template <typename Data>
   void GArray<Data>::append(const Data& data) 
   {
      assert(size_ <= capacity_);
      if (size_ == capacity_) {
         if (capacity_ == 0) {
            assert(data_ == 0); 
            Memory::allocate<Data>(data_, 64);
            capacity_ = 64;
         } else {
            assert(data_); 
            assert(capacity_ > 0); 
            Data* newPtr = 0;
            Memory::allocate<Data>(newPtr, 2*capacity_);
            if (size_ > 0) {
               for (int i = 0; i < size_; ++i) {
                  newPtr[i] = data_[i];
               }
            }
            Memory::deallocate<Data>(data_, capacity_);
            data_ = newPtr;
            capacity_ = 2*capacity_;
         }
      }
      // Append new element
      data_[size_] = data;
      ++size_;
      assert(size_ <= capacity_);
   }

   /*
   * Resize the array.
   */
   template <typename Data>
   void GArray<Data>::resize(int n) 
   {
      if (n < 0) {
         UTIL_THROW("Cannot resize to n < 0");
      }
      assert(capacity_ >= size_);
      if (n > size_) {
         if (n > capacity_) {
            int m = capacity_;
            if (m == 0) {
               m = n;
            } else {
               while (n > m) {
                  m *= 2;
               }
            }
            Data* newPtr = 0;
            Memory::allocate<Data>(newPtr, m);
            if (data_) {
               assert(capacity_ > 0);
               for (int i = 0; i < size_; ++i) {
                  newPtr[i] = data_[i];
               }
               Memory::deallocate<Data>(data_, capacity_);
            }
            data_ = newPtr;
            capacity_ = m;
         }
         // Initialize added elements, using placement new
         for (int i = size_; i < n;  ++i) {
            new(data_ + i) Data();
         }
      }
      size_ = n;
   }

   /*
   * Serialize a GArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void GArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      int size;
      if (Archive::is_saving()) {
         capacity = capacity_;
         size = size_;
      }
      ar & capacity;
      ar & size;
      if (Archive::is_loading()) {
         size_ = 0;
         reserve(capacity);
         size_ = size;
      }
      for (int i = 0; i < size_; ++i) {
         ar & data_[i];
      }
   }

   /*
   * Set an ArrayIterator to the beginning of this Array.
   */
   template <class Data>
   inline 
   void GArray<Data>::begin(ArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Set a ConstArrayIterator to the beginning of this Array.
   */
   template <class Data>
   inline 
   void GArray<Data>::begin(ConstArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline Data& GArray<Data>::operator[] (int i)
   {
      assert(i >= 0);
      assert(i < size_);
      return data_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline const Data& GArray<Data>::operator[] (int i) const
   {
      assert(i >= 0);
      assert(i < size_);
      return data_[i];
   }

   /*
   * Return physical capacity of array.
   */
   template <class Data>
   inline int GArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return logical size of this array (i.e., number of elements).
   */
   template <class Data>
   inline int GArray<Data>::size() const
   {  return size_; }

   /*
   * Is this array allocated?
   */
   template <class Data>
   inline bool GArray<Data>::isAllocated() const
   {  return (bool)data_; }

} 
#endif
