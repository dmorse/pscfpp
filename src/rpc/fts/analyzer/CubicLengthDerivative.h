#ifndef RPC_CUBIC_LENGTH_DERIVATIVE_H
#define RPC_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.h>  // base class template
#include <pscf/backends/CPT.h>                       // base class argument
#include <rp/fts/analyzer/AverageAnalyzer.h>       // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CubicLengthDerivative<1,CPT>;
      extern template class CubicLengthDerivative<2,CPT>;
      extern template class CubicLengthDerivative<3,CPT>;
   }
}
#endif
