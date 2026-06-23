/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Ramp.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rp/fts/ramp/Ramp.tpp>     // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor 
   template <int D>
   Ramp<D>::Ramp(Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : Rp::Ramp<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Ramp<1, Rpg::Types<1> >;
      template class Ramp<2, Rpg::Types<2> >;
      template class Ramp<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class Ramp<1>;
      template class Ramp<2>;
      template class Ramp<3>;
   }
}
