/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Perturbation.h"
#include <rpg/fts/simulator/Simulator.h>
#include <prdc/cuda/RField.h>

#include <rp/fts/perturbation/Perturbation.tpp> // base class implementation

namespace Pscf {
namespace Rpg {

   /*
   * Constructor.
   */
   template <int D>
   Perturbation<D>::Perturbation(Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : Rp::Perturbation<D, Types<D> > (simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Perturbation<1, Rpg::Types<1> >;
      template class Perturbation<2, Rpg::Types<2> >;
      template class Perturbation<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Perturbation<1>;
      template class Perturbation<2>;
      template class Perturbation<3>;
   }
}
