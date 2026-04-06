/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoxLengthDerivative.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/analyzer/BoxLengthDerivative.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   BoxLengthDerivative<D>::BoxLengthDerivative(Simulator<D>& simulator,
                                               System<D>& system)
    : Rp::BoxLengthDerivative<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class BoxLengthDerivative<1, Rpc::Types<1> >;
      template class BoxLengthDerivative<2, Rpc::Types<2> >;
      template class BoxLengthDerivative<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class BoxLengthDerivative<1>;
      template class BoxLengthDerivative<2>;
      template class BoxLengthDerivative<3>;
   }
}
