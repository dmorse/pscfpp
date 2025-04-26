#ifndef RPG_SIM_STATE_CU
#define RPG_SIM_STATE_CU

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.tpp"

namespace Pscf {
namespace Rpg {

   template struct SimState<1>;
   template struct SimState<2>;
   template struct SimState<3>;

}
}
#endif
