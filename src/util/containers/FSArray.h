#ifndef UTIL_FS_ARRAY_H
#define UTIL_FS_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ArrayIterator.h>
#include <util/containers/ConstArrayIterator.h>
#include <util/global.h>

namespace Util
{

   /**
   * A fixed capacity (static) contiguous array with a variable logical size.
   *
   * An FSArray < Data, Capacity > is a wrapper for a statically allocated 
   * C array containing Capacity objects of type Data. An FSArray has both
   * a Capacity that is set at compile time, which is the physical size of 
   * the underlying C array, and a logical size, which is the number of 
   * contiguous elements (from 0 to one less than its size) that contain 
   * valid data. The size is initialized to zero, and can only be increased 
   * only by the append() method, which adds a new element to the end of the 
   * array.  
   *
   * When compiled in debug mode (i.e., when NDEBUG is defined) the subcript
   * operator [] checks that the index is less than the logical size, and 
   * not merely less than the capacity.
   *
   * \ingroup Array_Module
   */
   template <typename Data, int Capacity>
   class FSArray 
   {

   public:

      /**
      * Constructor.
      */
      FSArray();

      /**
      * Copy constructor.
      *
      *\param other the FSArray to be copied.
      */
      FSArray(const FSArray<Data, Capacity>& other); 
   
      /**
      * Assignment, element by element.
      *
      * Capacity of LHS FSArray must be >= size of RHS FSArray.
      *
      * \param other the RHS FSArray 
      */
      FSArray<Data, Capacity>& operator=(const FSArray<Data, Capacity>& other);

      /**
      * Destructor.
      */
      virtual ~FSArray();

      /**
      * Return physical capacity of array.
      */
      int capacity() const;

      /**
      * Return logical size of this array (i.e., number of elements).
      */
      int size() const;

      /**
      * Set an ArrayIterator to the beginning of this container.
      *
      * \param iterator ArrayIterator, initialized on output. 
      */
      void begin(ArrayIterator<Data> &iterator);

      /**
      * Set a ConstArrayIterator to the beginning of this container.
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
      * Append data to the end of the array.
      *
      * \param data Data to add to end of array.
      */
      void append(const Data &data);

      /**
      * Set logical size to zero.
      */
      void clear();

      /**
      * Serialize to/from an archive.
      *
      * \param ar       archive 
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Packed size of FSArray in a MemoryArchive, in bytes.
      */
      int packedSize();

   protected:

      /// Array of Data elements.
      Data data_[Capacity];

      /// Logical size of array (number of elements used).
      int  size_;

   };

   /*
   * Constructor.
   */
   template <class Data, int Capacity>
   inline FSArray<Data, Capacity>::FSArray()
    : size_(0)
   {}

   /*
   * Copy constructor.
   *
   *\param other the FSArray to be copied.
   */
   template <class Data, int Capacity>
   FSArray<Data, Capacity>::FSArray(const FSArray<Data, Capacity>& other) 
   {
      size_  = other.size_;
      for (int i = 0; i < size_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment, element by element.
   *
   * Capacity of LHS FSArray must be >= size of RHS FSArray.
   *
   * \param other the RHS FSArray 
   */
   template <class Data, int Capacity>
   FSArray<Data, Capacity>& 
   FSArray<Data, Capacity>::operator=(const FSArray<Data, Capacity>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      // Copy elements
      size_ = other.size_;
      for (int i = 0; i < size_; ++i) {
         data_[i] = other[i];
      }
      return *this;
   }

   /*
   * Destructor.
   */
   template <class Data, int Capacity>
   FSArray<Data, Capacity>::~FSArray()
   {}

   /*
   * Return physical capacity of array.
   */
   template <class Data, int Capacity>
   int FSArray<Data, Capacity>::capacity() const
   { return Capacity; }

   /*
   * Return logical size of this array (i.e., number of elements).
   */
   template <class Data, int Capacity>
   int FSArray<Data, Capacity>::size() const
   { return size_; }

   /*
   * Set an ArrayIterator to the beginning of this Array.
   *
   * \param iterator ArrayIterator, initialized on output. 
   */
   template <class Data, int Capacity>
   void FSArray<Data, Capacity>::begin(ArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Set a ConstArrayIterator to the beginning of this Array.
   */
   template <class Data, int Capacity>
   void FSArray<Data, Capacity>::begin(ConstArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data, int Capacity>
   Data& FSArray<Data, Capacity>::operator[] (int i) 
   {
      assert(i < size_);
      assert(i >= 0);
      return data_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data, int Capacity>
   const Data& FSArray<Data, Capacity>::operator[] (int i) const
   {
      assert(i < size_);
      assert(i >= 0);
      return data_[i];
   }

   /*
   * Append data to the end of the array.
   */
   template <class Data, int Capacity>
   inline void FSArray<Data, Capacity>::append(const Data &data) 
   {
      if (size_ == Capacity) {
         UTIL_THROW("Attempt to add to full FSArray");
      }
      data_[size_] = data;
      ++size_;
   }

   /*
   * Set logical size to zero.
   */
   template <class Data, int Capacity>
   inline void FSArray<Data, Capacity>::clear() 
   {  size_ = 0; }

   /*
   * Serialize a FSArray to/from an Archive.
   */
   template <class Data, int Capacity>
   template <class Archive>
   inline void FSArray<Data, Capacity>::serialize(Archive& ar, 
                                          const unsigned int version)
   {
      ar & size_;
      if (size_ > Capacity) {
         UTIL_THROW("FSArray<Data, Capacity> with size > Capacity");
      }
      for (int i = 0; i < size_; ++i) {
         ar & data_[i];
      }
   }

   /*
   * Packed size of FSArray in a MemoryArchive, in bytes.
   */
   template <typename Data, int Capacity>
   inline int FSArray<Data, Capacity>::packedSize()
   {  return Capacity*sizeof(Data) + sizeof(int); }

} 
#endif
