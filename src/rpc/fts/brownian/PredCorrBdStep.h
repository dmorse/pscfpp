#ifndef RPC_PRED_CORR_BD_STEP_H
#define RPC_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/PredCorrBdStep.h> // base class template
#include <rpc/system/Types.h>               // base class argument 
#include <prdc/field/cpu/RField.h>          // base class member
#include <rpc/fts/brownian/BdStep.h>        // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::PredCorrBdStep<1, Rpc::Types<1> >;
      extern template class Rp::PredCorrBdStep<2, Rpc::Types<2> >;
      extern template class Rp::PredCorrBdStep<3, Rpc::Types<3> >;
   }
}
#endif
