/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StepLogger.h"
#include <rp/fts/analyzer/StepLogger.tpp>

namespace Pscf {
namespace Rpg {

   // Constructor.
   template <int D>
   StepLogger<D>::StepLogger(Rp::Simulator<D, Rpg::Types<D> >& simulator, Rp::System<D, Rpg::Types<D> >& system)
    : Rp::StepLogger< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class StepLogger< 1, Rpg::Types<1> >;
      template class StepLogger< 2, Rpg::Types<2> >;
      template class StepLogger< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class StepLogger<1>;
      template class StepLogger<2>;
      template class StepLogger<3>;
   }
}
