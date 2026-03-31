/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/environment/Environment.h>

#include <rp/scft/iterator/Iterator.tpp>

namespace Pscf {
namespace Rpc {

   /*
   * Default constructor.
   */
   template <int D>
   Iterator<D>::Iterator()
    : Rp::Iterator<D, System<D> >()
   {}

   /*
   * Constructor.
   */
   template <int D>
   Iterator<D>::Iterator(System<D>& system)
    : Rp::Iterator<D, System<D> >(system)
   {}

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Iterator<1, Rpc::System<1> >;
      template class Iterator<2, Rpc::System<2> >;
      template class Iterator<3, Rpc::System<3> >;
   }
   namespace Rpc {
      template class Iterator<1>;
      template class Iterator<2>;
      template class Iterator<3>;
   }
}
