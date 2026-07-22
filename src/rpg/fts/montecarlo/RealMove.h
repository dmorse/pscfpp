#ifndef RPG_REAL_MOVE_H
#define RPG_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/RealMove.h>    // base class template
#include <pscf/backends/CUT.h>              // template argument
#include <rpg/fts/montecarlo/McMove.h>     // indirect base class
#include <prdc/field/cuda/RField.h>               // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RealMove<1,CUT>;
      extern template class RealMove<2,CUT>;
      extern template class RealMove<3,CUT>;
   }
}
#endif
