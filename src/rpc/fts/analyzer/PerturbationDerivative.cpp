/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"
#include <rpc/fts/perturbation/Perturbation.h>
#include <rp/fts/analyzer/PerturbationDerivative.tpp>
#include <rpc/system/System.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <rpc/fts/simulator/Simulator.h>

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   PerturbationDerivative<D>::PerturbationDerivative(
                                     Simulator<D>& simulator,
                                     System<D>& system)
    : Rp::PerturbationDerivative< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationDerivative< 1, Rpc::Types<1> >;
      template class PerturbationDerivative< 2, Rpc::Types<2> >;
      template class PerturbationDerivative< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class PerturbationDerivative<1>;
      template class PerturbationDerivative<2>;
      template class PerturbationDerivative<3>;
   }
}
