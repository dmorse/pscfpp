/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/system/System.h>

#include <rp/fts/brownian/BdStep.tpp>

#if 0
namespace Pscf {
   namespace Rpc {
   
      // Constructor.
      template <int D>
      Rp::BdStep<D, Rpc::Types<D> >::BdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator)
       : Rp::BdStep<D, Types<D> >(simulator)
      {}
   
   }
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdStep< 1, Rpc::Types<1> >;
      template class BdStep< 2, Rpc::Types<2> >;
      template class BdStep< 3, Rpc::Types<3> >;
   }
}
