#ifndef PRDC_CUDA_FIELD_H
#define PRDC_CUDA_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Dynamic array on the GPU with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that is
   * allocated in device global memory.  All member functions may be called 
   * from the host, but the class thus does not offer access to individual 
   * elements via operator[]
   *
   * \ingroup Prdc_Cuda_Module
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
      Data* cField();

      /**
      * Return pointer to const to underlying C array.
      */
      const Data* cField() const;

      /**
      * Assignment operator.
      *
      * \param other Field<Data> on rhs of assignent (input)
      */
      virtual Field<Data>& operator = (const Field<Data>& other);

      /**
      * Copy constructor.
      * 
      * \param other Field<Data> to be copied (input)
      */
      Field(const Field& other);

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
   inline int Field<Data>::capacity() const
   {  return capacity_; }

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
   inline const Data* Field<Data>::cField() const
   {  return data_; }

   /*
   * Return true if the Field has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool Field<Data>::isAllocated() const
   {  return (bool)data_; }

}
}
}
#include "Field.tpp"
#endif
