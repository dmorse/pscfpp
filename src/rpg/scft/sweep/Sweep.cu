/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"

#include <rpg/system/System.h>
#include <rpg/scft/iterator/Iterator.h>
#include <rpg/scft/ScftThermo.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>

#include <rp/scft/sweep/Sweep.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Default constructor (for unit testing).
   template <int D>
   Sweep<D>::Sweep()
    : Rp::Sweep<D, Types<D> >()
   {}

   // Constructor, creates association with parent system.
   template <int D>
   Sweep<D>::Sweep(System<D> & sys)
    : Rp::Sweep<D, Types<D> >(sys)
   {}

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   template class SweepTmpl< Rpg::BasisFieldState<1> >;
   template class SweepTmpl< Rpg::BasisFieldState<2> >;
   template class SweepTmpl< Rpg::BasisFieldState<3> >;
   namespace Rp {
      template class Sweep<1, Rpg::Types<1> >;
      template class Sweep<2, Rpg::Types<2> >;
      template class Sweep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Sweep<1>;
      template class Sweep<2>;
      template class Sweep<3>;
   }
}
