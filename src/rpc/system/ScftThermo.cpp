/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftThermo.h"
#include <prdc/system/ScftReal.tpp>

namespace Pscf {
   namespace Prdc {
      template class ScftReal<1, Rpc::System<1> >;
      template class ScftReal<2, Rpc::System<2> >;
      template class ScftReal<3, Rpc::System<3> >;
   }
   namespace Rpc {
      template class ScftThermo<1>;
      template class ScftThermo<2>;
      template class ScftThermo<3>;
   }
}
