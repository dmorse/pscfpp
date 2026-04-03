/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>

#include <rp/scft/sweep/BasisFieldState.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   // Default constructor.
   template <int D>
   BasisFieldState<D>::BasisFieldState()
    : Rp::BasisFieldState<D, Types<D> > ()
   {}

   // Constructor - creates association with parent system.
   template <int D>
   BasisFieldState<D>::BasisFieldState(System<D>& system)
    : Rp::BasisFieldState<D, Types<D> > (system)
   {}


} // namespace Rpg
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BasisFieldState< 1, Rpg::Types<1> >;
      template class BasisFieldState< 2, Rpg::Types<2> >;
      template class BasisFieldState< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class BasisFieldState<1>;
      template class BasisFieldState<2>;
      template class BasisFieldState<3>;
   }
}
