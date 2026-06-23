/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepParameter.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/MixtureModifier.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>
#include <rpg/field/Domain.h>

#include <rp/scft/sweep/SweepParameter.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Default constructor.
   */
   template <int D>
   SweepParameter<D>::SweepParameter()
    : Rp::SweepParameter< D, Types<D> >()
   {}

   /*
   * Constructor, creates association with system.
   */
   template <int D>
   SweepParameter<D>::SweepParameter(Rp::System<D, Rpg::Types<D> >& system)
    : Rp::SweepParameter< D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1, Rpg::Types<1> >;
      template class SweepParameter<2, Rpg::Types<2> >;
      template class SweepParameter<3, Rpg::Types<3> >;
      
   }
   namespace Rpg {
      template class SweepParameter<1>;
      template class SweepParameter<2>;
      template class SweepParameter<3>;
   }
}
