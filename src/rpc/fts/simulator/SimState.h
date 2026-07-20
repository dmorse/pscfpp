#ifndef RPC_SIM_STATE_H
#define RPC_SIM_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/SimState.h>   // class template
#include <pscf/cpu/CppTp.h>            // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SimState<1, CppTp<1> >;
      extern template class SimState<2, CppTp<2> >;
      extern template class SimState<3, CppTp<3> >;
   }
}
#endif
