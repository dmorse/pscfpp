/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/montecarlo/McMoveManager.h>  // class header
#include <rpc/fts/montecarlo/McMoveFactory.h>
#include <rpc/fts/montecarlo/McSimulator.h>
#include <util/random/Random.h>
#include <util/global.h>

#include <rp/fts/montecarlo/McMoveManager.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveManager<1, CppTp<1> >;
      template class McMoveManager<2, CppTp<2> >;
      template class McMoveManager<3, CppTp<3> >;
   }
}
