/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CubicLengthDerivative.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <pscf/cpu/Reduce.h>

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   // Constructor.
   template <int D>
   CubicLengthDerivative<D>::CubicLengthDerivative(Simulator<D>& simulator,
                                               Rp::System<D, Rpc::Types<D> >& system)
    : Rp::CubicLengthDerivative<D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1, Rpc::Types<1> >;
      template class CubicLengthDerivative<2, Rpc::Types<2> >;
      template class CubicLengthDerivative<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class CubicLengthDerivative<1>;
      template class CubicLengthDerivative<2>;
      template class CubicLengthDerivative<3>;
   }
}
