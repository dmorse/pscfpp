/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "CPT.h"
#include <util/random/Random.h>

namespace Pscf {

   /*
   * Initialize backend.
   */
   void CPT::init()
   {}

   /*
   * Set thread count in backend.
   */
   void CPT::setThreadCount(int nThread)
   {}

   /*
   * Associate vector and scalar random number generators, if needed.
   *
   * GPU-enabled code would instead do nothing.
   */
   void CPT::linkVecRandom(VecRandom& vr, Random & sr)
   {  vr.associate(sr); }

   /*
   * Set vector RNG seed, if needed (do nothing for CPU code).
   *
   * GPU-enabled code would instead set the seed.
   */
   void CPT::seedVecRandom(VecRandom& vr, long seed)
   {}

}
