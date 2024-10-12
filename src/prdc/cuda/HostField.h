#ifndef PRDC_CUDA_HOST_FIELD_H
#define PRDC_CUDA_HOST_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Dynamic array of field values located on the host.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <typename Data>
   class HostField
   {

   public:

      /**
      * Default constructor.
      */
      HostField();

      /**
      * Copy constructor.
      * 
      * \param other HostField<Data> to be copied (input)
      */
      HostField(HostField<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~HostField();

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the HostField is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the HostField is not allocated.
      */
      virtual void deallocate();

      /**
      * Return true iff the HostField has been allocated.
      */
      bool isAllocated() const;

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
      * Assignment operator.
      *  
      * Performs a deep copy, by copying values of all elements.
      *
      * The RHS Field<Data> array must be allocated.  If this LHS 
      * HostField<D> is not allocated, memory will be allocated.  
      * If this LHS array is allocated, capacity values for LHS and 
      * RHS arrays must be equal. 
      *
      * \param other HostField<Data> on rhs of assignent (input)
      */
      virtual HostField<Data>& operator = (const HostField<Data>& other);

      /**
      * Assignment operator, assignment from Field<Data> device array.
      *
      * Performs a deep copy from Field<Data> RHS device array to this 
      * HostField<D> LHS host array, by copying underlying C array from
      * device memory to host memory.
      *
      * The RHS Field<Data> array must be allocated.  If this LHS 
      * HostField<D> is not allocated, allocate memory.  If this LHS
      * array is allocated, the capacity values for the LHS and RHS 
      * arrays must be equal.
      *
      * \param other Cpu::HostField<Data> on RHS of assignent (input)
      */
      virtual 
      HostField<Data>& operator = (const Cpu::HostField<Data>& other);

      /**
      * Return allocated size.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cArray();

      /**
      * Return pointer to const to underlying C array.
      */
      Data const * cArray() const;

   protected:

      /// Pointer to an array of Data elements.
      Data* data_;

      /// Allocated size of the data_ array.
      int capacity_;

   };

   /*
   * Get an element by reference (C-array subscripting)
   */
   template <typename Data>
   inline Data& HostField<Data>::operator[] (int i)
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
   inline Data const & HostField<Data>::operator[] (int i) const
   {
      assert(data_ != 0);
      assert(i >= 0 );
      assert(i < capacity_);
      return *(data_ + i);
   }

   /*
   * Return allocated size.
   */
   template <typename Data>
   inline int HostField<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* HostField<Data>::cArray()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline 
   Data const * HostField<Data>::cArray() const
   {  return data_; }

   /*
   * Return true if the HostField has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool HostField<Data>::isAllocated() const
   {  return (bool) data_; }

}
}
}
#include "HostField.tpp"
#endif
