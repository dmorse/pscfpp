#ifndef PSCF_CUT_H
#define PSCF_CUT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/BackendId.h>
#include <pscf/backend/cuda/cudaTypes.h>

// Forward declarations
namespace Util {
   class Random;
}
namespace Pscf {
   template <typename T> class VecRandom;
}

namespace Pscf {

   using namespace Util;

   /**
   * Backend identifier class for the GPU-enabled CUDA backend.
   *
   * \ingroup Pscf_Backend_Cuda_Module
   */
   class CUT {

   public:

      // Type aliases

      using Real = cudaReal;
      using Complex = cudaComplex;

      // Static members

      /**
      * Identifier for computational backend.
      */
      static const BackendId id = BackendId::Cpp;

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

   };

} // namespace Pscf
#endif
