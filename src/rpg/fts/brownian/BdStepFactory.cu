/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
// Subclasses of BdStep
#include <rp/fts/brownian/ExplicitBdStep.h>
#include <rp/fts/brownian/PredCorrBdStep.h>
#include <rp/fts/brownian/LMBdStep.h>
#endif

#include <rp/fts/brownian/BdStepFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdStepFactory<1,CUT>;
      template class BdStepFactory<2,CUT>;
      template class BdStepFactory<3,CUT>;
   }
}
