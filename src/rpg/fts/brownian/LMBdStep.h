#ifndef RPG_LM_BD_STEP_H
#define RPG_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/LMBdStep.h>     // base class template
#include <rpg/system/Types.h>             // base class template argument 
#include <rpg/fts/brownian/BdStep.h>      // indirect base class
#include <prdc/field/cuda/RField.h>       // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::LMBdStep<1, Rpg::Types<1> >;
      extern template class Rp::LMBdStep<2, Rpg::Types<2> >;
      extern template class Rp::LMBdStep<3, Rpg::Types<3> >;
   }
}
#endif
