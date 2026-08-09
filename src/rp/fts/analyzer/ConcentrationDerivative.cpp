/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/fts/analyzer/ConcentrationDerivative.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationDerivative<1,CPT>;
      template class ConcentrationDerivative<2,CPT>;
      template class ConcentrationDerivative<3,CPT>;
   }
}
