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
   static void Types<D>::init()
   {  ThreadArray::init(); }

   /*
   * Set thread count in backend.
   */
   template <int D>
   static void Types<D>::setThreadCount(int nThread)
   {
      ThreadArray::setThreadsPerBlock(nThread);
      ThreadMesh::setThreadsPerBlock(nThread);
   }

}
}
