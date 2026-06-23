/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/WFields.h>

#include <rp/scft/sweep/BasisFieldState.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   // Default constructor.
   template <int D>
   BasisFieldState<D>::BasisFieldState()
    : Rp::BasisFieldState<D, Types<D> > ()
   {}

   // Constructor - creates association with parent system.
   template <int D>
   BasisFieldState<D>::BasisFieldState(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::BasisFieldState<D, Types<D> > (system)
   {}


} // namespace Rpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BasisFieldState< 1, Rpc::Types<1> >;
      template class BasisFieldState< 2, Rpc::Types<2> >;
      template class BasisFieldState< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class BasisFieldState<1>;
      template class BasisFieldState<2>;
      template class BasisFieldState<3>;
   }
}
