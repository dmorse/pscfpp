#ifndef RPG_SIM_STATE_H
#define RPG_SIM_STATE_H

/*
*PSCF - Polymer Self-Consistent Field
*
*Copyright 2015 - 2025, The Regents of the University of Minnesota
*Distributed under the terms of the GNU General Public License.
*/

#include<rp/fts/simulator/SimState.h>  // class template
#include<rpg/system/Types.h>           // template argument

//Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template struct SimState<1, Rpg::Types<1> >;
      extern template struct SimState<2, Rpg::Types<2> >;
      extern template struct SimState<3, Rpg::Types<3> >;
   }
}
#endif
