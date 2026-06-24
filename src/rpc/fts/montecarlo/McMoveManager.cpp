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

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   McMoveManager<D>::McMoveManager(Rp::McSimulator<D, Rpc::Types<D> >& simulator,
                                   Rp::System<D, Rpc::Types<D> >& system)
    : Rp::McMoveManager<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McMoveManager<1, Rpc::Types<1> >;
      template class McMoveManager<2, Rpc::Types<2> >;
      template class McMoveManager<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class McMoveManager<1>;
      template class McMoveManager<2>;
      template class McMoveManager<3>;
   }
}
