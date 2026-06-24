#ifndef RPG_SWEEP_H
#define RPG_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/Sweep.h>             // base class template
#include <rpg/system/Types.h>                // base class argument
#include <rpg/scft/sweep/BasisFieldState.h>  // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template class SweepTmpl< Rp::BasisFieldState<1, Rpg::Types<1> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<2, Rpg::Types<2> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<3, Rpg::Types<3> > >;
   namespace Rp {
      extern template class Sweep<1, Rpg::Types<1> >;
      extern template class Sweep<2, Rpg::Types<2> >;
      extern template class Sweep<3, Rpg::Types<3> >;
   }
}
#endif
