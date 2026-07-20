/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StepLogger.h"
#include <rp/fts/analyzer/StepLogger.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class StepLogger< 1, CudaTp<1> >;
      template class StepLogger< 2, CudaTp<2> >;
      template class StepLogger< 3, CudaTp<3> >;
   }
}
