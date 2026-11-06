/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/fts/analyzer/ConcentrationDerivative.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationDerivative<1,CPT>;
      template class ConcentrationDerivative<2,CPT>;
      template class ConcentrationDerivative<3,CPT>;
   }
}
