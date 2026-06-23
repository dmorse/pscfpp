/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampParameter.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/system/System.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/MixtureModifier.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>
#include <rpg/field/Domain.h>

#include <rp/fts/ramp/RampParameter.tpp>    // base class implementation

namespace Pscf {
namespace Rpg {

   // Default constructor.
   template <int D>
   RampParameter<D>::RampParameter()
    : Rp::RampParameter<D, Types<D> >()
   {}

   // Constructor - creates association with a Simulator.
   template <int D>
   RampParameter<D>::RampParameter(Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : Rp::RampParameter<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampParameter<1, Rpg::Types<1> >;
      template class RampParameter<2, Rpg::Types<2> >;
      template class RampParameter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class RampParameter<1>;
      template class RampParameter<2>;
      template class RampParameter<3>;
   }
}
