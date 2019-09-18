#ifndef PSSP_GPU_RING_BUFFER_H
#define PSSP_GPU_RING_BUFFER_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * Class for storing history of previous values in an array.
   *
   * Data is stored in a circular buffer, in which, once the
   * array is full, the newest data value overwrites the oldest.
   *
   * \ingroup Array_Module
   */
   template <class Data>
   class RingBuffer
   {
   
   public: 
   
      /**
      * Contructor.
      */
      RingBuffer();
   
      /**
      * Copy contructor.
      *
      * \param other object to be copied.
      */
      RingBuffer(const RingBuffer<Data>& other);
   
      /**
      * Assignment.
      *
      * \param other object to be assigned.
      */
      RingBuffer& operator = (const RingBuffer<Data>& other);
   
      /**
      * Destructor.
      */
      virtual ~RingBuffer();
   
      /**
      * Allocate a new empty buffer.
      *
      * Throw an Exception if this RingBuffer has already been
      * allocated - a RingBuffer can only be allocated once.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Set previously allocated buffer to empty state.
      */
      void clear();
  
      /**
      * Add a new value to the buffer.
      *
      * \param value new value to be added.
      */
      void append(Data value);
   
      /**
      * Return number of values currently in the buffer.
      */
      int size() const;
  
      /**
      * Return the capacity of the buffer.
      */
      int capacity() const;
  
      /** 
      * Return true if the RingBuffer has been allocated, false otherwise.
      */
      bool isAllocated() const;
   
      /** 
      * Return true if full (if size == capacity), false otherwise.
      */
      bool isFull();
   
      /**
      * Retrieve a previous value, a specified number of time steps ago.
      *
      * \param offset number of steps back in time (offset=0 is current value).
      */
      Data& operator [] (int offset) const;

      /**
      * Serialize a RingBuffer to/from an Archive.
      *
      * \param ar      archive 
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:
  
      /// Underlying C-array of Data elements.
      Data* data_;     

      /// Allocated size of data_ array (zero before allocation).
      int   capacity_; 

      /// Number of values in buffer, must be <= capacity_.
      int   size_;     

      /// Integer array index of the most recently added value.
      int   last_;     
  
   };

   /*
   * Default contructor
   */
   template <class Data>
   RingBuffer<Data>::RingBuffer()
    : data_(0),
      capacity_(0),
      size_(0),
      last_(0)
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RingBuffer to be copied.
   */
   template <class Data>
   RingBuffer<Data>::RingBuffer(const RingBuffer<Data>& other) 
    : data_(0),
      capacity_(0),
      size_(0),
      last_(0)
   {
      if (other.capacity_ > 0) {
         assert(other.data_ != 0);
         Memory::allocate<Data>(data_, other.capacity_);
         capacity_ = other.capacity_;
         size_ = other.size_;
         last_ = other.last_;
         for (int i = 0; i < capacity_; ++i) {
            data_[i] = other.data_[i];
         }
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * Assignment is allowed only for allocated RingBuffer objects of equal 
   * capacity.  This operator will allocate this RingBuffer if it was not 
   * allocated previously.
   */
   template <class Data>
   RingBuffer<Data>& 
   RingBuffer<Data>::operator = (const RingBuffer<Data>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      if (other.isAllocated()) {

         assert(other.capacity_ > 0);

         if (!isAllocated()) {

            Memory::allocate<Data>(data_, other.capacity_);
            capacity_ = other.capacity_;

         } else if (capacity_ != other.capacity_) {

            UTIL_THROW("Cannot assign RingBuffers of unequal capacity");

         }

         size_ = other.size_;
         last_ = other.last_;

         for (int i = 0; i < capacity_; ++i) {
            data_[i] = other.data_[i];
         }

      } else if (isAllocated()) {
         UTIL_THROW("Cannot assign unallocated array to allocated array");
      }

      return *this;
   }

   /// Destructor
   template <class Data>
   RingBuffer<Data>::~RingBuffer()
   {
      if (data_) {
         Memory::deallocate<Data>(data_, capacity_);
      }
   }

   /// Allocate a new array containing capacity elements.
   template <class Data>
   void  RingBuffer<Data>::allocate(int capacity) 
   {
      if (data_ == 0) {
         Memory::allocate<Data>(data_, capacity);
         capacity_ = capacity;
      } else {
         UTIL_THROW("Error: Attempt to re-allocate a RingBuffer");
      }
      last_ = capacity_ - 1; 
      size_ = 0;        // No values in buffer
   }

   /*
   * Set buffer to empty state, by resetting counters.
   */
   template <class Data>
   void RingBuffer<Data>::clear()
   {
      last_  = capacity_ - 1; 
      size_  = 0;        
   }
 
   /* 
   * Append a new value to the end of the buffer
   */
   template <class Data>
   void RingBuffer<Data>::append(Data data) 
   {
      assert(last_ < capacity_);
      ++last_;
      if (last_ == capacity_) {
         last_ = 0; // wrap around
      };
      data_[last_] = data;
   
      // If buffer is not yet full, update size_
      if (size_ < capacity_) {
         ++size_;
      };

   }
   
   /*
   * Return number of values in buffer.
   */
   template <class Data>
   int RingBuffer<Data>::size() const
   { return size_; }
  
   /*
   * Return the capacity of the buffer.
   */
   template <class Data>
   int RingBuffer<Data>::capacity() const
   { return capacity_; }
 
   /* 
   * Return true if the RingBuffer has been allocated, false otherwise.
   */
   template <class Data>
   bool RingBuffer<Data>::isAllocated() const
   {  return (data_ != 0); }
  
   /* 
   * Return true if the RingBuffer is full.
   */
   template <class Data>
   bool RingBuffer<Data>::isFull()
   {  return (size_ == capacity_); }
  
   /*
   * Retrive a previous value, index backwards from 0 (current)
   */
   template <class Data>
   Data& RingBuffer<Data>::operator [] (int offset) const
   {
      if ( offset >= size_ )  {
         UTIL_THROW("offset >= size_");
      }
      int i;
      i = last_ - offset;
      if (i < 0) {
         i = i + capacity_;
      }
      return data_[i];
   }

   /*
   * Serialize a RingBuffer to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void RingBuffer<Data>::serialize(Archive& ar, const unsigned int version)
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
            } else {
               capacity_ = 0;
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent RingBuffer capacities");
            }
         }
      }
      ar & size_;
      ar & last_;
      for (int i = 0; i < capacity_; ++i) {
         ar & data_[i];
      }
   }

}
#endif
