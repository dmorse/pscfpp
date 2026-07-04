#ifndef RPC_CUBIC_LENGTH_DERIVATIVE_H
#define RPC_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.h>
#include <rpc/system/Types.h>
#include <rpc/fts/analyzer/AverageAnalyzer.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CubicLengthDerivative<1, Rpc::Types<1> >;
      extern template class CubicLengthDerivative<2, Rpc::Types<2> >;
      extern template class CubicLengthDerivative<3, Rpc::Types<3> >;
   }
}
#endif
