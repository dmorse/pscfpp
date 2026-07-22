#ifndef RPG_SWEEP_H
#define RPG_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/Sweep.h>             // base class template
#include <pscf/backends/CUT.h>                // base class argument
#include <rpg/scft/sweep/BasisFieldState.h>  // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template class SweepTmpl< Rp::BasisFieldState<1,CUT> >;
   extern template class SweepTmpl< Rp::BasisFieldState<2,CUT> >;
   extern template class SweepTmpl< Rp::BasisFieldState<3,CUT> >;
   namespace Rp {
      extern template class Sweep<1,CUT>;
      extern template class Sweep<2,CUT>;
      extern template class Sweep<3,CUT>;
   }
}
#endif
