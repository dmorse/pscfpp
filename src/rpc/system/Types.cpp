/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "Types.h"
#include <pscf/cpu/CpuVecRandom.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   /*
   * Initialize backend.
   */
   template <int D>
   void Types<D>::init()
   {}

   /*
   * Set thread count in backend.
   */
   template <int D>
   void Types<D>::setThreadCount(int nThread)
   {}

   /*
   * Associate vector and scalar random number generators, if needed.
   *
   * GPU-enabled code would instead do nothing.
   */
   template <int D>
   void Types<D>::linkVecRandom(VecRandom& vr, Random & sr)
   {  vr.associate(sr); }

   /*
   * Set vector RNG seed, if needed (do nothing for CPU code).
   *
   * GPU-enabled code would instead set the seed.
   */
   template <int D>
   void Types<D>::seedVecRandom(VecRandom& vr, long seed)
   {}

   // Explicit instantiation definitions
   template class Types<1>;
   template class Types<2>;
   template class Types<3>;

}
}
