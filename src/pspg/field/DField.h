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
   * Dynamic array on the GPU with alligned data.
   *
   * This class wraps an aligned C array with elements of type Data on the 
   * device. All member functions may be called from the host. As a result,
   * the class does not offer access to individual elements via operator[]
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
      * Allocate the underlying C array on the device.
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

      /**
      * Assignment operator.
      *
      * \param other DField<Data> on rhs of assignent (input)
      */
      virtual DField<Data>& operator = (const DField<Data>& other);

      /**
      * Copy constructor.
      * 
      * \param other DField<Data> to be copied (input)
      */
      DField(const DField& other);

   protected:

      /// Pointer to an array of Data elements on the device / GPU.
      Data* data_;

      /// Allocated size of the data_ array.
      int capacity_;

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

}
}
#include "DField.tpp"
#endif
