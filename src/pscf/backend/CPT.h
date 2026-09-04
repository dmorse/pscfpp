#ifndef PSCF_CPT_H
#define PSCF_CPT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/BackendId.h>
#include <pscf/cpu/FftwDRArray.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <fftw3.h>

// Forward declarations
namespace Util {
   class Random;
}

namespace Pscf {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * Type class for the C++ serial CPU backend.
   */
   class CPT
   {

   public:

      // Aliases for associated types

      using Real = double;
      using Complex = fftw_complex;

      using VecRandom = CpuVecRandom;

      template <typename T> using DevArray = FftwDRArray<T>;
      template <typename T> using LocArray = FftwDRArray<T>;
      using RDevArray = DevArray<Real>;
      using RLocArray = LocArray<Real>;

      // Static members

      /**
      * Identifier for computational backend.
      */
      static const Prdc::BackendId id = BackendId::Cpp;

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

      /**
      * Link vector and scalar random number generators, if needed.
      *
      * The implementation for CPU code associates the two generators.
      * In CPU code, the vector RNG simply calls the scalar RNG.
      *
      * \param vr vector RNG
      * \param sr scalar RNG
      */
      static
      void linkVecRandom(VecRandom& vr, Random & sr);

      /**
      * Set the vector random number generator seed, if needed.
      *
      * The implementation for CPU code does nothing, because the vector
      * RNG is simply a facade for the scalar RNG.
      *
      * \param vr  vector RNG to be initialized
      * \param seed  random seed value
      */
      static
      void seedVecRandom(VecRandom& vr, long seed);

   };

} // namespace Pscf
#endif
