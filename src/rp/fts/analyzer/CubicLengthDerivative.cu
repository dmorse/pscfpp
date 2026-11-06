/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1,CUT>;
      template class CubicLengthDerivative<2,CUT>;
      template class CubicLengthDerivative<3,CUT>;
   }
}
