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
   * GPU specializatin allocate host array if not done previously. 
   *
   * \param in  input array to be copied
   * \param out  output array into which data is copied
   */
   template <typename T>
   void setupHostArray(HostDArray<T> & hostArray, 
		       DeviceArray<T> const & deviceArray)
   {
      if (!hostArray.isAllocated()) {
         int n = deviceArray.capacity();
         hostArray.allocate(n);
      }
   }

   /**
   * Copy data from device to host.
   *
   * This specialization for GPU hardware actually copies data.
   *
   * \param in  input array to be copied
   * \param out  output array into which data is copied
   */
   template <typename T>
   void sendToHost(HostDArray<T>& out, DeviceArray<T> const & in)
   {  out = in; }

   /**
   * Copy data from host to device.
   *
   * This specialization for GPU hardware actually copies data.
   *
   * \param in  input array to be copied
   * \param out  output array into which data is copied
   */
   template <typename T>
   void sendToDevice(DeviceArray<T> & out, HostDArray<T> const & in)
   {  out = in; }

   /**
   * Release host array for re-use.
   *
   * This specialization for GPU hardware does nothing. The host array
   * will be be de-allocated when it goes out of scope.
   *
   * \param array  array to be released
   */
   template <typename T>
   void releaseHostArray(HostDArray<T> & array)
   {}

}
#endif
