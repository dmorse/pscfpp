/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"

#include <rpc/system/System.h>
#include <rpc/scft/iterator/Iterator.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/scft/sweep/Sweep.tpp>

namespace Pscf {
namespace Rpc {

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

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   template class SweepTmpl< Rpc::BasisFieldState<1> >;
   template class SweepTmpl< Rpc::BasisFieldState<2> >;
   template class SweepTmpl< Rpc::BasisFieldState<3> >;
   namespace Rp {
      template class Sweep<1, Rpc::Types<1> >;
      template class Sweep<2, Rpc::Types<2> >;
      template class Sweep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class Sweep<1>;
      template class Sweep<2>;
      template class Sweep<3>;
   }
}
