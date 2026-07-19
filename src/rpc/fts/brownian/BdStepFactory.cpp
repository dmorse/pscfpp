/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStepFactory.h"

// Subclasses of BdStep
#include <rpc/fts/brownian/ExplicitBdStep.h>
#include <rpc/fts/brownian/PredCorrBdStep.h>
#include <rpc/fts/brownian/LMBdStep.h>
//#include <rpc/fts/brownian/BdSimulator.h>

#include <rp/fts/brownian/BdStepFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdStepFactory<1, Cpp<1> >;
      template class BdStepFactory<2, Cpp<2> >;
      template class BdStepFactory<3, Cpp<3> >;
   }
}
