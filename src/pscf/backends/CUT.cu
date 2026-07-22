/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "CUT.h"
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/ThreadMesh.h>
#include <util/random/Random.h>

namespace Pscf {

   /*
   * Initialize backend.
   */
   void CUT::init()
   {  ThreadArray::init(); }

   /*
   * Set thread count in backend.
   */
   void CUT::setThreadCount(int nThread)
   {
      ThreadArray::setThreadsPerBlock(nThread);
      ThreadMesh::setThreadsPerBlock(nThread);
   }

   /*
   * Link vector and scalar random number generators, if needed.
   *
   * The Cuda GPU implementation does nothing, because host and device 
   * RNGs are are independent for this code. The CPU implementation would
   * link the two RNGs.
   */
   void CUT::linkVecRandom(VecRandom& vr, Random & sr)
   {}

   /*
   * Set vector RNG seed, if needed (do nothing for CPU code).
   *
   * The Cuda GPU implementation sets the seed for the GPU vector RNG.
   * The CPU implementation would do nothing, because the RNGs are linked
   * in this case.
   */
   void CUT::seedVecRandom(VecRandom& vr, long seed)
   {  vr.setSeed(seed); }

}
