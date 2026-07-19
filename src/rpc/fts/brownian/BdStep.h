#ifndef RPC_BD_STEP_H
#define RPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>
#include <pscf/cpu/Cpp.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdStep< 1, Cpp<1> >;
      extern template class BdStep< 2, Cpp<2> >;
      extern template class BdStep< 3, Cpp<3> >;
   }
}
#endif
