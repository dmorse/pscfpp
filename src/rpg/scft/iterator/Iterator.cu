/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/environment/Environment.h>

#include <rp/scft/iterator/Iterator.tpp>

namespace Pscf {
namespace Rpg {

   /*
   * Default constructor.
   */
   template <int D>
   Iterator<D>::Iterator()
    : Rp::Iterator<D, Types<D> >()
   {}

   /*
   * Constructor.
   */
   template <int D>
   Iterator<D>::Iterator(System<D>& system)
    : Rp::Iterator<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Iterator<1, Rpg::Types<1> >;
      template class Iterator<2, Rpg::Types<2> >;
      template class Iterator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Iterator<1>;
      template class Iterator<2>;
      template class Iterator<3>;
   }
}
