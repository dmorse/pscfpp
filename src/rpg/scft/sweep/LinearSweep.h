#ifndef RPG_LINEAR_SWEEP_H
#define RPG_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/LinearSweep.h>      // direct base class template
#include <pscf/backends/CUT.h>               // base class argument
#include <rpg/scft/sweep/Sweep.h>           // indirect base class
#include <rpg/scft/sweep/SweepParameter.h>  // indirect base member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LinearSweep<1,CUT>;
      extern template class LinearSweep<2,CUT>;
      extern template class LinearSweep<3,CUT>;
   }
}
#endif
