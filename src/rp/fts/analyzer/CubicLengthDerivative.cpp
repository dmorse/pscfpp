/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>
#include <pscf/backend/cpp/CPT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1,CPT>;
      template class CubicLengthDerivative<2,CPT>;
      template class CubicLengthDerivative<3,CPT>;
   }
}
