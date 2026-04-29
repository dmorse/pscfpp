/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiDerivative.h"                    // header
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/fts/analyzer/ChiDerivative.tpp>  // base class implementation

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   ChiDerivative<D>::ChiDerivative(Simulator<D>& simulator,
                                   System<D>& system)
    : Rp::ChiDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ChiDerivative< 1, Rpc::Types<1> >;
      template class ChiDerivative< 2, Rpc::Types<2> >;
      template class ChiDerivative< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ChiDerivative<1>;
      template class ChiDerivative<2>;
      template class ChiDerivative<3>;
   }
}
