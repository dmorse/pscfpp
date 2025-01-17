#ifndef PSCF_HOST_D_ARRAY_H
#define PSCF_HOST_D_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   // Forward declaration of analogous container for data on the device.
   template <typename Data> class DeviceArray;

   /**
   * Template for dynamic array stored in host CPU memory.
   *
   * This class is provided as a convenience to allow the use of 
   * assigment (=) operators to copy data between corresponding 
   * containers that store array data in device vs. host memory. A 
   * HostDArray<Data> stores data in a dynamically allocated array in 
   * host CPU memory, whereas a DeviceArray<Data> stores analogous 
   * data in global GPU device memory. Each of these classes defines  
   * an assigment operation that allows assignment from the other, 
   * which silently copies the underlying arrays between device and 
   * host memory.
   *
   * The underlying array is allocated using cudaMallocHost.
   *
   * \ingroup Pscf_Cuda_Module
   */
   template <typename Data>
   class HostDArray
   {

   public:

      /**
      * Default constructor.
      */
      HostDArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      HostDArray(int capacity);

      /**
      * Copy constructor.
      * 
      * \param other HostDArray<Data> to be copied (input)
      */
      HostDArray(HostDArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~HostDArray();

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the HostDArray is already allocated.
      *
      * \param capacity number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the HostDArray is not allocated
      */
      virtual void deallocate();

      /**
      * Assignment operator, assignment from another HostDArray<Data>.
      *  
      * Performs a deep copy, by copying all elements of the underlying
      * C array from host memory to host memory.
      *
      * The RHS HostDArray<Data> object must be allocated.  If this LHS 
      * HostDArray<D> is not allocated, the correct size block of memory 
      * will be allocated. Otherwise, if this LHS object is allocated
      * on entry, capacity values for LHS and RHS objects must be equal. 
      *
      * \throw Exception if other HostDArray<Data> is not allocated
      * Exceptions are thrown if the RHS array is not allocated or if
      * both arrays are allocated but have unequal capacities.
      *
      * \param other HostDArray<Data> on RHS of assignment (input)
      */
      virtual
      HostDArray<Data>& operator = (const HostDArray<Data>& other);

      /**
      * Assignment operator, assign from DeviceArray<Data> device array.
      *
      * Performs a deep copy from a RHS DeviceArray<Data> device array 
      * to this LHS HostDArray<D> host array, by copying the underlying 
      * C array from device memory to host memory.
      *
      * The RHS DeviceArray<Data> object must be allocated.  If this LHS 
      * HostDArray<D> is not allocated, the correct size block of memory 
      * will be allocated. Otherwise, if this LHS object is allocated on
      * entry, capacity values for LHS and RHS objects must be equal. 
      *
      * Exceptions are thrown if the RHS array is not allocated or if
      * both arrays are allocated but have unequal capacities.
      *
      * \param other DeviceArray<Data> array on RHS of assignment (input)
      */
      virtual 
      HostDArray<Data>& operator = (const DeviceArray<Data>& other);

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
      Data* cArray();

      /**
      * Return pointer to const to the underlying C array on the host.
      */
      Data const * cArray() const;

      /**
      * Return true iff the HostDArray has been allocated.
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
   inline Data& HostDArray<Data>::operator[] (int i)
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
   inline Data const & HostDArray<Data>::operator[] (int i) const
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
   inline int HostDArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* HostDArray<Data>::cArray()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline 
   Data const * HostDArray<Data>::cArray() const
   {  return data_; }

   /*
   * Return true if the HostDArray has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool HostDArray<Data>::isAllocated() const
   {  return (bool) data_; }

   #ifndef PSCF_HOST_D_ARRAY_TPP
   extern template class HostDArray<cudaReal>;
   extern template class HostDArray<cudaComplex>;
   extern template class HostDArray<int>;
   extern template class HostDArray<bool>;
   #endif

}
#endif
