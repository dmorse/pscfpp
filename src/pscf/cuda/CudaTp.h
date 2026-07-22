#ifndef PSCF_CUDA_TP_H
#define PSCF_CUDA_TP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/BackendId.h>
#include <pscf/cuda/cudaTypes.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>

// Forward declarations
namespace Util {
   class Random;
}
namespace Pscf {
   class CudaVecRandom;
}

namespace Pscf {

   // Namespace that may be used implicitly
   using namespace Util;

   /**
   * Type class for the Cuda GPU backend.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D>
   class CudaTp {

   public:

      // Type aliases

      using Real = cudaReal;
      using Complex = cudaComplex;

      template <typename T> using DevArray = DeviceArray<T>;
      template <typename T> using LocArray = HostDArray<T>;

      using RDevArray = DevArray<Real>;
      using RLocArray = LocArray<Real>;

      using VecRandom = CudaVecRandom;

      // Static members

      /**
      * Identifier for computational backend.
      */
      static const Prdc::BackendId id = Prdc::BackendId::Cpp;

      /**
      * Initialize backend thread array.
      */
      static void init();

      /**
      * Set the number of threads per block.
      *
      * \param nThread  number of threads
      */
      static void setThreadCount(int nThread);

      /**
      * Link vector and scalar random number generators, if needed.
      *
      * The Cuda GPU implementation does nothing, because the host and
      * device RNGs are independent.
      *
      * \param vr  vector RNG
      * \param sr  scalar RNG
      */
      static
      void linkVecRandom(VecRandom& vr, Random & sr);

      /**
      * Set the vector random number generator seed, if needed.
      *
      * The Cuda GPU implementation sets a seed for the device RNG.
      *
      * \param vr  vector RNG to be initialized
      * \param seed  random seed value
      */
      static
      void seedVecRandom(VecRandom& vr, long seed);

   };

   // Explicit instantiation declarations
   extern template class CudaTp<1>;
   extern template class CudaTp<2>;
   extern template class CudaTp<3>;

} // namespace Pscf
#endif
