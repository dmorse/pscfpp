#ifndef PRDC_CUDA_SEND_H
#define PRDC_CUDA_SEND_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/DeviceArray.h> 
#include <pscf/cuda/HostDArray.h> 

namespace Pscf {

   using namespace Util;

   /**
   * Setup host array for use.
   *
   * GPU specialization allocates host array if not done previously. 
   *
   * \param hostArray host array that should be allocated if necessary
   * \param deviceArray  device array, which must be allocated on entry
   */
   template <typename T>
   void setupHostArray(HostDArray<T> & hostArray, 
		       DeviceArray<T> const & deviceArray)
   {
      UTIL_CHECK(deviceArray.isAllocated());
      const int n = deviceArray.capacity();
      if (!hostArray.isAllocated()) {
         hostArray.allocate(n);
      }
      UTIL_CHECK(hostArray.capacity() == n);
   }

   /**
   * Copy data from device to host.
   *
   * This specialization for GPU hardware actually copies data from 
   * GPU device memory to CPU host memory.
   *
   * \param in  input device array to be copied
   * \param out  output host array into which data is copied
   */
   template <typename T>
   void sendToHost(HostDArray<T>& out, DeviceArray<T> const & in)
   {  out = in; }

   /**
   * Copy data from host to device.
   *
   * This specialization for GPU hardware actually copies data from
   * CPU host memory to GPU device memory.
   *
   * \param in  input host array to be copied
   * \param out  output device array into which data is copied
   */
   template <typename T>
   void sendToDevice(DeviceArray<T> & out, HostDArray<T> const & in)
   {  out = in; }

   /**
   * Release host array for re-use.
   *
   * This specialization for GPU hardware does nothing. The host array
   * will be be de-allocated when it goes out of scope. This pattern
   * allows use of a host array that is a class member that is only
   * allocated once.
   *
   * \param array  host array to be released
   */
   template <typename T>
   void releaseHostArray(HostDArray<T> & array)
   {}

}
#endif
