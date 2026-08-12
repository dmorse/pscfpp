#ifndef CPC_SYSTEM_CPP
#define CPC_SYSTEM_CPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"              // header
#include <cpc/solvers/Mixture.h>
#include <cpc/field/Domain.h>
#include <prdc/field/cpu/WaveList.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/CField.h>
#include <pscf/interaction/Interaction.h>

#include <cp/system/System.tpp>  // base class template implementation 

namespace Pscf {

   namespace Cp {
      // Explicit instantiation of base class
      template class System< 1, Cpc::Types<1> >;
      template class System< 2, Cpc::Types<2> >;
      template class System< 3, Cpc::Types<3> >;
   }

   namespace Cpc {
      // Explicit instantiation
      template class System<1>;
      template class System<2>;
      template class System<3>;
   }

}
#endif
