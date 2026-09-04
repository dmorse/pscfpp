/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/Reduce.h>
#include <pscf/backend/cpp/CPT.h>

#include <rp/fts/analyzer/CubicLengthDerivative.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class CubicLengthDerivative<1,CPT>;
      template class CubicLengthDerivative<2,CPT>;
      template class CubicLengthDerivative<3,CPT>;
   }
}
