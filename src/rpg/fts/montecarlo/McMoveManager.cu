/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/fts/montecarlo/McMoveManager.h>  // class header
#include <rpg/fts/montecarlo/McMoveFactory.h>
#include <rpg/fts/montecarlo/McSimulator.h>
#include <util/random/Random.h>
#include <util/global.h>

#include <rp/fts/montecarlo/McMoveManager.tpp> // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   McMoveManager<D>::McMoveManager(McSimulator<D>& simulator,
                                   Rp::System<D, Rpg::Types<D> >& system)
    : Rp::McMoveManager<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveManager<1, Rpg::Types<1> >;
      template class McMoveManager<2, Rpg::Types<2> >;
      template class McMoveManager<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class McMoveManager<1>;
      template class McMoveManager<2>;
      template class McMoveManager<3>;
   }
}
