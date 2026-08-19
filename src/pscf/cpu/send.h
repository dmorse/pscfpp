#ifndef PSCF_CPU_SEND_H
#define PSCF_CPU_SEND_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/FftwDRArray.h> 
#include <util/global.h> 

namespace Pscf {

   using namespace Util;

   /*
   * These functions are designed to allow a pair of FftwDArray containers 
   * to be used in template code in a manner analogous to the use of a
   * device and host array in code that stores memory in the separate 
   * memory of a GPU device. Function template specializations designed 
   * for use with a GPU actually copy memory between the device and host. 
   * The specializations defined here, which are designed to use only a
   * CPU, instead create a temporary association between the two arrays 
   * without copying any data. In this pattern, the "device array" is the 
   * array that owns the underlying C array, while the "host" array is a 
   * shallow copy that is associated with the device array. 
   */

   /**
   * Setup host/user array for use.
   *
   * CPU specialization associates host array with a device array.
   *
   * \param hostArray  host array (data user)
   * \param deviceArray  device array (data owner)
   */
   template <typename T>
   void setupHostArray(FftwDRArray<T> & hostArray, 
                       FftwDRArray<T> & deviceArray)
   {
      UTIL_CHECK(deviceArray.isOwner());
      UTIL_CHECK(!hostArray.isOwner());
      if (!hostArray.isAssociated()) {
         const int begin = 0;
         const int capacity = deviceArray.capacity();
         hostArray.associate(deviceArray, begin, capacity);
      }
      UTIL_CHECK(hostArray.isAssociated());
   }

   /**
   * Copy data from device array to host array.
   *
   * This specialization for CPU hardware does nothing, since the two
   * arrays must be associated on entry to this function. An Exception is 
   * thrown if the input is not a data owner, or if the output host array
   * is not already associated with the input.
   *
   * \param in  input device array
   * \param out  output host array 
   */
   template <typename T>
   void sendToHost(FftwDRArray<T>& out, FftwDRArray<T> const & in)
   {
      UTIL_CHECK(in.isOwner());
      UTIL_CHECK(out.isAssociated());
   }

   /**
   * Copy data from host array to device array.
   *
   * This specialization for CPU hardware does nothing, since the two
   * arrays must be associated on entry to this function. An Exception is 
   * thrown if the output device array is not a data owner, or if the 
   * input host array is not already associated with the output.
   *
   * \param in  input host array
   * \param out  output device array
   */
   template <typename T>
   void sendToDevice(FftwDRArray<T> & out, FftwDRArray<T> const & in)
   {
      UTIL_CHECK(out.isOwner());
      UTIL_CHECK(in.isAssociated());
   }

   /**
   * Release host array for re-use.
   *
   * This specialization for CPU hardware dissociates the user host array
   * from an associated data owner.
   *
   * \param array  host array to be released
   */
   template <typename T>
   void releaseHostArray(FftwDRArray<T> & array)
   {
      UTIL_CHECK(!array.isOwner());
      UTIL_CHECK(array.isAssociated());
      array.dissociate();
   }

}
#endif
