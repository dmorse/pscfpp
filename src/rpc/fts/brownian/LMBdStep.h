#ifndef RPC_LM_BD_STEP_H
#define RPC_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/LMBdStep.h>      // base class template
#include <pscf/cpu/CppTp.h>              // base class argument 
#include <rpc/fts/brownian/BdStep.h>       // indirect base class
#include <prdc/field/cpu/RField.h>         // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LMBdStep<1, CppTp<1> >;
      extern template class LMBdStep<2, CppTp<2> >;
      extern template class LMBdStep<3, CppTp<3> >;
   }
}
#endif
