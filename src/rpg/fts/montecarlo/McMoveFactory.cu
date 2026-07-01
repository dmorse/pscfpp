/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveFactory.h"

// Subclasses of McMove
#include <rpg/fts/montecarlo/RealMove.h>
#include <rpg/fts/montecarlo/ForceBiasMove.h>
#include <rpg/fts/montecarlo/BdMove.h>
#include <rpg/fts/montecarlo/ShiftMove.h>

#include <rp/fts/montecarlo/McMoveFactory.tpp> // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveFactory<1, Rpg::Types<1> >;
      template class McMoveFactory<2, Rpg::Types<2> >;
      template class McMoveFactory<3, Rpg::Types<3> >;
   }
}
