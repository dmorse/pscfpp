/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/


#include "CUT.h"
#include "VecRandom.h"
#include <pscf/backend/cuda/ThreadArray.h>
#include <pscf/backend/cuda/ThreadMesh.h>
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

}
