#ifndef PSPG_DFIELD_H
#define PSPG_DFIELD_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /**
   * Dynamic array with aligned data, for use with cufftw library/device code.
   * This class does not offer memory access via operator[]
   *
   * \ingroup Pspg_Field_Module
   */
   template <typename Data>
   class DField
   {

   public:

      /**
      * Default constructor.
      */
      DField();

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DField();

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
      * Return pointer to underlying C array.
      */
      Data* cDField();

      /**
      * Return pointer to const to underlying C array.
      */
      const Data* cDField() const;


      //Removing this. Child class has this function
      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      /*
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);
      */

   protected:

      /// Pointer to an array of Data elements.
      Data* data_;

      /// Allocated size of the data_ array.
      int capacity_;

   private:

      /**
      * Copy constructor (private and not implemented to prohibit).
      */
      DField(const DField& other);

      /**
      * Assignment operator (private and non implemented to prohibit).
      */
      DField& operator = (const DField& other);

   };

   /*
   * Return allocated size.
   */
   template <typename Data>
   inline int DField<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* DField<Data>::cDField()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline const Data* DField<Data>::cDField() const
   {  return data_; }

   /*
   * Return true if the Field has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool DField<Data>::isAllocated() const
   {  return (bool)data_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   /*template <typename Data>
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
   }*/

}
}
#include "DField.tpp"
#endif
