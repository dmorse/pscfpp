#ifndef RPG_SWEEP_PARAMETER_H
#define RPG_SWEEP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/SweepParameter.h>
#include <rpg/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SweepParameter<1, Rpg::Types<1> >;
      extern template class SweepParameter<2, Rpg::Types<2> >;
      extern template class SweepParameter<3, Rpg::Types<3> >;
   }
}
#endif
