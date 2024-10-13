#ifndef PRDC_CUDA_HOST_FIELD_H
#define PRDC_CUDA_HOST_FIELD_H

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

   // Forward declaration of analogous container for data on the device.
   template <typename Data> class Field;

   /**
   * Dynamic array of field values in host CPU memory.
   *
   * This class is provided as a convenience to allow the use of assigment (=)
   * operators to copy data between corresponding containers that store field
   * data in device vs. host memory. A HostField<Data> stores data in a 
   * dynamically allocated array in host memory, whereas a Field<Data> stores 
   * analogous data in global device memory. Each of these classes defines an 
   * assigment operation that allows assignment from the other, and that 
   * silently copies the underlying arrays between device and host memory.
   *
   * Memory is allocated using cudaMallocHost.
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
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      HostField(int capacity);

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
      * \param capacity number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the HostField is not allocated
      */
      virtual void deallocate();

      /**
      * Assignment operator, assignment from another HostField<Data>.
      *  
      * Performs a deep copy, by copying all elements of the underlying
      * C array from host memory to host memory.
      *
      * The RHS HostField<Data> object must be allocated.  If this LHS 
      * HostField<D> is not allocated, the correct size block of memory 
      * will be allocated. Otherwise, if this LHS object is allocated
      * on entry, capacity values for LHS and RHS objects must be equal. 
      *
      * \throw Exception if other HostField<Data> is not allocated
      * Exceptions are thrown if the RHS array is not allocated or if
      * both arrays are allocated but have unequal capacities.
      *
      * \param other HostField<Data> on RHS of assignment (input)
      */
      virtual
      HostField<Data>& operator = (const HostField<Data>& other);

      /**
      * Assignment operator, assignment from Field<Data> device array.
      *
      * Performs a deep copy from a RHS Field<Data> device array to this 
      * LHS HostField<D> host array, by copying the underlying C array 
      * from device memory to host memory.
      *
      * The RHS Field<Data> object must be allocated.  If this LHS 
      * HostField<D> is not allocated, the correct size block of memory 
      * will be allocated. Otherwise, if this LHS object is allocated on
      * entry, capacity values for LHS and RHS objects must be equal. 
      *
      * Exceptions are thrown if the RHS array is not allocated or if
      * both arrays are allocated but have unequal capacities.
      *
      * \param other Field<Data> device array on RHS of assignment (input)
      */
      virtual 
      HostField<Data>& operator = (const Field<Data>& other);

      /**
      * Get one element by non-const reference.
      *
      * Mimic C-array subscripting.
      *
      * \param  i array index
      * \return non-const reference to element i
      */
      Data & operator[] (int i);

      /**
      * Get one element by const reference.
      *
      * Mimics C-array subscripting.
      *
      * \param i array index
      * \return const reference to element i
      */
      Data const & operator[] (int i) const;

      /**
      * Return allocated array capacity.
      *
      * \return Number of elements allocated in the array
      */
      int capacity() const;

      /**
      * Return pointer to the underlying C array on the host.
      */
      Data* cField();

      /**
      * Return pointer to const to the underlying C array on the host.
      */
      Data const * cField() const;

      /**
      * Return true iff the HostField has been allocated.
      */
      bool isAllocated() const;

   private:

      /// Pointer to a C array of Data elements, allocated on host.
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
   * Return allocated array capacity.
   */
   template <typename Data>
   inline int HostField<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* HostField<Data>::cField()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline 
   Data const * HostField<Data>::cField() const
   {  return data_; }

   /*
   * Return true if the HostField has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool HostField<Data>::isAllocated() const
   {  return (bool) data_; }

   #ifndef PRDC_CUDA_HOST_FIELD_TPP
   extern template class HostField<cudaReal>;
   extern template class HostField<cudaComplex>;
   #endif

}
}
}
#endif
