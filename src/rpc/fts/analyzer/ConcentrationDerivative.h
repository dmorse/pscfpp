#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rp/fts/analyzer/ConcentrationDerivative.h>
#include <rpc/system/Types.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationDerivative<1, Rpc::Types<1> >;
      extern template class ConcentrationDerivative<2, Rpc::Types<2> >;
      extern template class ConcentrationDerivative<3, Rpc::Types<3> >;
   }
}
#endif
