#ifndef RPG_RAMP_PARAMETER_H
#define RPG_RAMP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/RampParameter.h>   // direct base class template
#include <rpg/system/Types.h>            // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RampParameter<1, Rpg::Types<1> >;
      extern template class RampParameter<2, Rpg::Types<2> >;
      extern template class RampParameter<3, Rpg::Types<3> >;
   }
}
#endif
