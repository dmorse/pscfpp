#ifndef PSCF_CPT_H
#define PSCF_CPT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/BackendId.h>
#include <fftw3.h>

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
   * Backend identifier class for the C++ serial CPU backend.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   class CPT
   {

   public:

      // Aliases for associated types

      using Real = double;
      using Complex = fftw_complex;

      // Static members

      /**
      * Identifier for computational backend.
      */
      static const BackendId id = BackendId::Cpp;

      /**
      * Initialize backend.
      */
      static void init();

      /**
      * Set the number of threads.
      *
      * \param nThread  number of threads
      */
      static void setThreadCount(int nThread);

   };

} // namespace Pscf
#endif
