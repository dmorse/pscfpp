/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "Types.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/ThreadMesh.h>

namespace Pscf {
namespace Rpg {

   /*
   * Initialize backend.
   */
   template <int D>
   void Types<D>::init()
   {  ThreadArray::init(); }

   /*
   * Set thread count in backend.
   */
   template <int D>
   void Types<D>::setThreadCount(int nThread)
   {
      ThreadArray::setThreadsPerBlock(nThread);
      ThreadMesh::setThreadsPerBlock(nThread);
   }

   /*
   * Initialize the random number generator.
   */
   template <int D>
   void Types<D>::initVecRandom(VecRandom& vr, Random & r, long seed)
   {  vr.setSeed(seed); }

   /*
   * Link vector and scalar random number generators, if needed.
   *
   * GPU implementation does nothing, because host and device RNGs are
   * are independent in this case.
   */
   template <int D>
   void Types<D>::linkVecRandom(VecRandom& vr, Random & sr)
   {}

   /*
   * Set vector RNG seed, if needed (do nothing for CPU code).
   *
   * GPU implementation sets the seed for the GPU vector RNG.
   */
   template <int D>
   void Types<D>::seedVecRandom(VecRandom& vr, long seed)
   {  vr.setSeed(seed); }

   // Explicit instantiation definitions
   template class Types<1>;
   template class Types<2>;
   template class Types<3>;

}
}
