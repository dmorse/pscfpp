/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CubicLengthDerivative.h"

#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/CUT.h>

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1,CUT>;
      template class CubicLengthDerivative<2,CUT>;
      template class CubicLengthDerivative<3,CUT>;
   }
}
