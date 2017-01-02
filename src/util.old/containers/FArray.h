#ifndef UTIL_F_ARRAY_H
#define UTIL_F_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ArrayIterator.h>
#include <util/containers/ConstArrayIterator.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>   
#include <util/mpi/MpiStructBuilder.h>   
#endif 

namespace Util
{

   /**
   * A fixed size (static) contiguous array template.
   *
   * An FArray is a simple wraper for a fixed size C Array, with a
   * capacity that is fixed at compile time. As in a C Array, or a
   * DArray container, all of the elements are accessible. Unlike
   * an FSArray, an FArray does not have logical size that is 
   * distinct from its physical capacity.
   *
   * When bounds checking is on (i.e., when NDEBUG is not defined),
   * the operator [] checks that the index is non-negative and less
   * than the Capacity.
   *
   * Advice: Use an FArray if you know exactly how many elements will 
   * be needed at compile time. Use an FSArray when you need a small 
   * statically allocated array for which the maximum capacity needed
   * is known at compile time, but the logical size may be less than 
   * the capacity. Use a DArray if you need a large, dynamically 
   * allocated array that must be allocated after instantiation.
   *
   * \ingroup Array_Module
   */
   template <typename Data, int Capacity>
   class FArray 
   {

   public:

      /**
      * Constructor.
      */
      FArray();

      /**
      * Copy constructor.
      *
      *\param other the FArray to be copied.
      */
      FArray(const FArray<Data, Capacity>& other);
   
      /**
      * Assignment, element by element.
      *
      * Capacity of LHS FArray must be >= size of RHS FArray.
      *
      * \param other the RHS FArray 
      */
      FArray<Data, Capacity>& operator=(const FArray<Data, Capacity>& other);

      // Default destructor is okay.

      /**
      * Return number of elements in this FArray.
      */
      int size() const;

      /**
      * Return number of elements in this FArray.
      */
      int capacity() const;

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
      * Return pointer to underlying C array.
      */
      Data* cArray();

      /**
      * Return pointer to const to underlying C array.
      */
      const Data* cArray() const;

      /**
      * Serialize a FArray to/from an Archive.
      *
      * \param ar        archive 
      * \param version   archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Return packed size in a MemoryArchive, in bytes.
      */
      int packedSize();

      #ifdef UTIL_MPI
      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();
      #endif

   private:

      /// Array of Data elements.
      Data data_[Capacity];

   };

   /*
   * Constructor.
   */
   template <typename Data, int Capacity>
   FArray<Data, Capacity>::FArray()
   {}

   /*
   * Copy constructor.
   *
   *\param other the FArray to be copied.
   */
   template <typename Data, int Capacity>
   FArray<Data, Capacity>::FArray(const FArray<Data, Capacity>& other) 
   {
      for (int i = 0; i < Capacity; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment, element by element.
   *
   * Capacity of LHS FArray must be >= size of RHS FArray.
   *
   * \param other the RHS FArray 
   */
   template <typename Data, int Capacity>
   FArray<Data, Capacity>& 
   FArray<Data, Capacity>::operator=(const FArray<Data, Capacity>& other) 
   {

      // Check for self assignment
      if (this == &other) return *this;

      // Copy elements
      for (int i = 0; i < Capacity; ++i) {
         data_[i] = other[i];
      }
      return *this;
   }

   // Default constructor is okay.

   /*
   * Return number of elements in this FArray.
   */
   template <typename Data, int Capacity>
   inline int FArray<Data, Capacity>::size() const
   { return Capacity; }

   /*
   * Return number of elements in this FArray.
   */
   template <typename Data, int Capacity>
   inline int FArray<Data, Capacity>::capacity() const
   { return Capacity; }

   /*
   * Set an ArrayIterator to the beginning of this Array.
   */
   template <typename Data, int Capacity>
   inline void FArray<Data, Capacity>::begin(ArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + Capacity);
   }

   /*
   * Set a ConstArrayIterator to the beginning of this Array.
   */
   template <typename Data, int Capacity>
   inline void 
   FArray<Data, Capacity>::begin(ConstArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + Capacity);
   }

   /*
   * Mimic C array subscripting.
   */
   template <typename Data, int Capacity>
   inline Data& FArray<Data, Capacity>::operator[] (int i)
   {
      assert(i < Capacity);
      assert(i >= 0);
      return data_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <typename Data, int Capacity>
   inline const Data& FArray<Data, Capacity>::operator[] (int i) const
   {
      assert(i < Capacity);
      assert(i >= 0 );
      return data_[i];
   }

   /*
   * Return pointer to underlying C array.
   */
   template <typename Data, int Capacity>
   inline Data* FArray<Data, Capacity>::cArray()
   {  return data_; }

   /*
   * Return pointer to const to underlying C array.
   */
   template <typename Data, int Capacity>
   const Data* FArray<Data, Capacity>::cArray() const
   {  return data_; }

   /*
   * Serialize a FArray to/from an Archive.
   */
   template <class Data, int Capacity>
   template <class Archive>
   void FArray<Data, Capacity>::serialize(Archive& ar, 
                                          const unsigned int version)
   {
      for (int i = 0; i < Capacity; ++i) {
         ar & data_[i];
      }
   }

   /**
   * Packed size of FArray in a MemoryArchive, in bytes.
   */
   template <typename Data, int Capacity>
   int FArray<Data, Capacity>::packedSize()
   {  return Capacity*sizeof(Data); }

   #ifdef UTIL_MPI
   /*
   * Commit associated MPI Datatype.
   */
   template <typename Data, int Capacity>
   void FArray<Data, Capacity>::commitMpiType() 
   {
      if (!MpiTraits< FArray<Data, Capacity > >::hasType) {
         MpiStructBuilder       builder;
         FArray<Data, Capacity> object;
   
         builder.setBase(&object);
         for (int i = 0; i < Capacity; ++i) {
            builder.addMember(&(object.data_[i]), MpiTraits<Data>::type);
         }
         builder.commit(MpiTraits< FArray<Data, Capacity> >::type);
         MpiTraits< FArray<Data, Capacity> >::hasType = true;
      }
   }
   #endif 
} 
#endif
