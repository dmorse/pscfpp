/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveFactory.h"

// Subclasses of McMove
#include <rpc/fts/montecarlo/RealMove.h>
#include <rpc/fts/montecarlo/ForceBiasMove.h>
#include <rpc/fts/montecarlo/BdMove.h>
#include <rpc/fts/montecarlo/ShiftMove.h>

#include <rp/fts/montecarlo/McMoveFactory.tpp> // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveFactory<1, CppTp<1> >;
      template class McMoveFactory<2, CppTp<2> >;
      template class McMoveFactory<3, CppTp<3> >;
   }
}
