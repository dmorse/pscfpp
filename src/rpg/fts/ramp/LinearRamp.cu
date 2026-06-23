/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearRamp.h"
#include <rpg/fts/simulator/Simulator.h>
#include <rp/fts/ramp/LinearRamp.tpp>      // base class implementation

namespace Pscf {
namespace Rpg {

   // Constructor
   template <int D>
   LinearRamp<D>::LinearRamp(Rp::Simulator<D, Rpg::Types<D> >& simulator)
    : Rp::LinearRamp<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class LinearRamp<1, Rpg::Types<1> >;
      template class LinearRamp<2, Rpg::Types<2> >;
      template class LinearRamp<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class LinearRamp<1>;
      template class LinearRamp<2>;
      template class LinearRamp<3>;
   }
}
