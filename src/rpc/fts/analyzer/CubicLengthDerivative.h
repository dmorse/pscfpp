#ifndef RPC_CUBIC_LENGTH_DERIVATIVE_H
#define RPC_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.h>  // base class template
#include <pscf/cpu/Cpp.h>                       // base class argument
#include <rpc/fts/analyzer/AverageAnalyzer.h>       // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CubicLengthDerivative<1, Cpp<1> >;
      extern template class CubicLengthDerivative<2, Cpp<2> >;
      extern template class CubicLengthDerivative<3, Cpp<3> >;
   }
}
#endif
