/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "Types.h"

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

   // Explicit instantiation definitions
   template class Types<1>;
   template class Types<2>;
   template class Types<3>;

}
}
