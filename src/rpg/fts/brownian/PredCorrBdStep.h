#ifndef RPG_PRED_CORR_BD_STEP_H
#define RPG_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/PredCorrBdStep.h> // base class template
#include <pscf/backends/CUT.h>               // base class template argument 
#include <prdc/field/cuda/RField.h>               // base class member
#include <rpg/fts/brownian/BdStep.h>        // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::PredCorrBdStep<1,CUT>;
      extern template class Rp::PredCorrBdStep<2,CUT>;
      extern template class Rp::PredCorrBdStep<3,CUT>;
   }
}
#endif
