/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"            // class header

#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <pscf/cuda/Reduce.h>
#include <prdc/cuda/RField.h>
#include <pscf/interaction/Interaction.h>

#include <rp/scft/ScftThermo.tpp>  // base class template implementation

namespace Pscf {
   namespace Rpg {

      // Constructor
      template <int D>
      ScftThermo<D>::ScftThermo(System<D> const & system)
       : Rp::ScftThermo<D, System<D> >(system)
      {};

   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1, Rpg::System<1> >;
      template class ScftThermo<2, Rpg::System<2> >;
      template class ScftThermo<3, Rpg::System<3> >;
   }
   namespace Rpg {
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;

   }
}
