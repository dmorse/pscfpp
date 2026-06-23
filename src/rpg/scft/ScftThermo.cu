/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"                     // class header
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/cuda/Reduce.h>

#include <rp/scft/ScftThermo.tpp>           // base class implementation

#if 0
namespace Pscf {
   namespace Rpg {

      /*
      * Constructor.
      */
      template <int D>
      Rp::ScftThermo<D, Rpg::Types<D> >::ScftThermo(Rp::System<D, Rpg::Types<D> > const & system)
       : Rp::ScftThermo<D, Types<D> >(system)
      {};

   }
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ScftThermo<1, Rpg::Types<1> >;
      template class ScftThermo<2, Rpg::Types<2> >;
      template class ScftThermo<3, Rpg::Types<3> >;
   }
   #if 0
   namespace Rpg {
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;
   }
   #endif
}
