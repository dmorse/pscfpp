#ifndef RPC_CHI_DERIVATIVE_H
#define RPC_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                // indirect base class
#include <rp/fts/analyzer/BinaryChiDerivative.h>  // base class template
#include <rpc/system/Types.h>               // base template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryChiDerivative<1, Rpc::Types<1> >;
      extern template class BinaryChiDerivative<2, Rpc::Types<2> >;
      extern template class BinaryChiDerivative<3, Rpc::Types<3> >;
   }
}
#endif
