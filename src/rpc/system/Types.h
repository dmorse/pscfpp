#ifndef RPC_TYPES_H
#define RPC_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>
#include <prdc/field/BackendId.h>

// Forward declarations
namespace Util {
   template <typename T> class DArray;
   template <typename T> class DRArray;
   class Random;
}
namespace Pscf {
   class CpuVecRandom;
   namespace Prdc {
      namespace Cpu {
         template <typename T> class FftwDRArray;
         template <int D> class RField;
         template <int D> class RFieldDft;
         template <int D> class FFT;
         template <int D> class RFieldComparison;
         template <int D> class RFieldDftComparison;
         template <int D> class WaveList;
      }
   }
}

namespace Pscf {
namespace Rpc {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * Aliases for types used by the C++ Cpu backend.
   */
   template <int D>
   class Types 
   {

   public:

      // Aliases for associated types

      using Real = double;
      using Complex = fftw_complex;

      using VecRandom = CpuVecRandom;

      template <typename T> using DevArray = Prdc::Cpu::FftwDRArray<T>;
      template <typename T> using LocArray = Prdc::Cpu::FftwDRArray<T>;
      using RDevArray = DevArray<Real>;
      using RLocArray = LocArray<Real>;

      using RField = Prdc::Cpu::RField<D>;
      using RFieldDft = Prdc::Cpu::RFieldDft<D>;
      using FFT = Prdc::Cpu::FFT<D>;
      using RFieldComparison = Prdc::Cpu::RFieldComparison<D>;
      using RFieldDftComparison = Prdc::Cpu::RFieldDftComparison<D>;
      using WaveList = Prdc::Cpu::WaveList<D>;

      // Static members

      /**
      * Identifier for computational backend.
      */
      static const Prdc::BackendId backendId = BackendId::Cpp;

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

   // Explicit instantiation declarations
   extern template class Types<1>;
   extern template class Types<2>;
   extern template class Types<3>;

} 
} // namespace Pscf
#endif
