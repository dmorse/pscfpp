#ifndef PSPC_FIELD_H
#define PSPC_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * Dynamic array with aligned data, for use with FFTW library.
   *
   * \ingroup Pspc_Field_Module
   */
   template <typename Data>
   class Field
   {

   public:

      /**
      * Default constructor.
      */
      Field();

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~Field();

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the Field is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the Field is not allocated.
      */
      void deallocate();

      /**
      * Return true if the Field has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Return allocated size.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Get an element by non-const reference.
      *
      * Mimic C-array subscripting.
      *
      * \param  i array index
      * \return non-const reference to element i
      */
      Data & operator[] (int i);

      /**
      * Get an element by const reference.
      *
      * Mimics C-array subscripting.
      *
      * \param i array index
      * \return const reference to element i
      */
      Data const & operator[] (int i) const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cField();

      /**
      * Return pointer to const to underlying C array.
      */
      Data const * cField() const;

      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   protected:

      /// Pointer to an array of Data elements.
      Data* data_;

      /// Allocated size of the data_ array.
      int capacity_;

   private:

      /**
      * Copy constructor (private and not implemented to prohibit).
      */
      Field(Field const & other);

      /**
      * Assignment operator (private and non implemented to prohibit).
      */
      Field& operator = (Field const & other);

   };

   /*
   * Return allocated size.
   */
   template <typename Data>
   inline int Field<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get an element by reference (C-array subscripting)
   */
   template <typename Data>
   inline Data& Field<Data>::operator[] (int i)
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity_);
      return *(data_ + i);
   }

   /*
   * Get an element by const reference (C-array subscripting)
   */
   template <typename Data>
   inline Data const & Field<Data>::operator[] (int i) const
   {
      assert(data_ != 0);
      assert(i >= 0 );
      assert(i < capacity_);
      return *(data_ + i);
   }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* Field<Data>::cField()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline 
   Data const * Field<Data>::cField() const
   {  return data_; }

   /*
   * Return true if the Field has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool Field<Data>::isAllocated() const
   {  return (bool)data_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <typename Data>
   template <class Archive>
   void Field<Data>::serialize(Archive& ar, const unsigned int version)
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
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent Field capacities");
            }
         }
      }
      if (isAllocated()) {
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

}
}
#include "Field.tpp"
#endif
