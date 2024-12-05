#ifndef PSCF_DEVICE_D_ARRAY_H
#define PSCF_DEVICE_D_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GpuTypes.h"
#include <util/global.h>
#include <cufft.h>

namespace Pscf {

   using namespace Util;

   // Forward declaration of analogous container for data on host.
   template <typename Data> class HostDArray;

   /**
   * Dynamic array on the GPU device with aligned data.
   *
   * This class wraps an aligned C array with elements of type Data that 
   * is allocated in GPU device global memory.  All member functions may 
   * be called from the CPU host, but this class does not offer access 
   * to individual elements via the subscript operator, operator[].
   *
   * \ingroup Pscf_Cuda_Module
   */
   template <typename Data>
   class DeviceDArray
   {

   public:

      /**
      * Default constructor.
      */
      DeviceDArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      * 
      * \param capacity number of elements to allocate 
      */
      DeviceDArray(int capacity);

      /**
      * Copy constructor.
      * 
      * \param other DeviceDArray<Data> to be copied (input)
      */
      DeviceDArray(DeviceDArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DeviceDArray();

      /**
      * Allocate the underlying C array on the device.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the array is not allocated.
      */
      void deallocate();

      /**
      * Assignment operator, assign from another DeviceDArray<Data> array.
      *  
      * Performs a deep copy, by copying values of all elements from 
      * device memory to device memory.
      *
      * This function will allocate memory if this (LHS) array is not 
      * allocated.  If this is allocated, it must have the same 
      * dimensions as the RHS DeviceDArray<Data>.
      *
      * \param other DeviceDArray<Data> on rhs of assignent (input)
      */
      virtual 
      DeviceDArray<Data>& operator = (const DeviceDArray<Data>& other);

      /**
      * Assignment operator, assignment from HostDArray<Data> host array.
      *
      * Performs a deep copy from a RHS HostDArray<Data> host array to 
      * this LHS DeviceDArray<Data> device array, by copying underlying 
      * C array from host memory to device memory.
      *
      * This function will allocate memory if this (LHS) 
      * DeviceDArray<Data> is not allocated.  If this is allocated, it 
      * must have the same dimensions as the RHS HostDArray<Data>.
      *
      * \param other HostDArray<Data> on RHS of assignent (input)
      */
      virtual 
      DeviceDArray<Data>& operator = (const HostDArray<Data>& other);

      /**
      * Return allocated capacity.
      *
      * \return Number of elements allocated in array
      */
      int capacity() const;

      /**
      * Return true if the array has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Return pointer to underlying C array.
      */
      Data* cArray();

      /**
      * Return pointer to const to underlying C array.
      */
      Data const * cArray() const;

   protected:

      /// Pointer to a C array of Data elements on the GPU device.
      Data* data_;

      /// Allocated size (capacity) of the data_ array.
      int capacity_;

   };

   /*
   * Return allocated capacity.
   */
   template <typename Data>
   inline int DeviceDArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get a pointer to the underlying C array.
   */
   template <typename Data>
   inline Data* DeviceDArray<Data>::cArray()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   template <typename Data>
   inline const Data* DeviceDArray<Data>::cArray() const
   {  return data_; }

   /*
   * Return true if the array has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool DeviceDArray<Data>::isAllocated() const
   {  return (bool)data_; }

   #ifndef PSCF_DEVICE_D_ARRAY_TPP
   extern template class DeviceDArray<cudaReal>;
   extern template class DeviceDArray<cudaComplex>;
   #endif

} // namespace Pscf
#endif
