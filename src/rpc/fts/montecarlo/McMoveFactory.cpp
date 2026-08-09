/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveFactory.h"

// Subclasses of McMove
#include <rp/fts/montecarlo/RealMove.h>
#include <rp/fts/montecarlo/ForceBiasMove.h>
#include <rp/fts/montecarlo/BdMove.h>
#include <rp/fts/montecarlo/ShiftMove.h>

#include <rp/fts/montecarlo/McMoveFactory.tpp> // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveFactory<1,CPT>;
      template class McMoveFactory<2,CPT>;
      template class McMoveFactory<3,CPT>;
   }
}
