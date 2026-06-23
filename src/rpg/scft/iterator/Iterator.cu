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

#if 0
namespace Pscf {
namespace Rpg {

   /*
   * Default constructor.
   */
   template <int D>
   Rp::Iterator<D, Rpg::Types<D> >::Iterator()
    : Rp::Iterator<D, Types<D> >()
   {}

   /*
   * Constructor.
   */
   template <int D>
   Rp::Iterator<D, Rpg::Types<D> >::Iterator(Rp::System<D, Rpg::Types<D> >& system)
    : Rp::Iterator<D, Types<D> >(system)
   {}

}
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Iterator<1, Rpg::Types<1> >;
      template class Iterator<2, Rpg::Types<2> >;
      template class Iterator<3, Rpg::Types<3> >;
   }
}
