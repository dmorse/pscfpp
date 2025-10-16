/*
* PSCF - Mixture Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"
#include <prdc/cpu/FFT.h>

namespace Pscf {
   namespace Prdc { 
      template class MixturePrdc<1, Rpc::Polymer<1>, Rpc::Solvent<1>, Rpc::Types<1> >;
      template class MixturePrdc<2, Rpc::Polymer<2>, Rpc::Solvent<2>, Rpc::Types<2> >;
      template class MixturePrdc<3, Rpc::Polymer<3>, Rpc::Solvent<3>, Rpc::Types<3> >;
   }
   namespace Rpc { 
      template class Mixture<1>;
      template class Mixture<2>;
      template class Mixture<3>;
   }
}
