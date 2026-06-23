#ifndef RPG_SIMULATOR_H
#define RPG_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/Simulator.h>    // base class template
#include <rpg/system/Types.h>              // template argument
#include <rpg/fts/simulator/SimState.h>    // member
#include <prdc/cuda/RField.h>              // member (template arg)

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Simulator<1, Rpg::Types<1> >;
      extern template class Simulator<2, Rpg::Types<2> >;
      extern template class Simulator<3, Rpg::Types<3> >;
   }
}
#endif
