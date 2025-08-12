#ifndef RPC_SIM_STATE_CPP
#define RPC_SIM_STATE_CPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.tpp"

namespace Pscf {
namespace Rpc {

   template struct SimState<1>;
   template struct SimState<2>;
   template struct SimState<3>;

}
}
#endif
