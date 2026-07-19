#ifndef RPC_CHI_DERIVATIVE_H
#define RPC_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryChiDerivative.h>  // base class template
#include <pscf/cpu/Cpp.h>                     // base template argument
#include <rpc/fts/analyzer/AverageAnalyzer.h>     // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryChiDerivative<1, Cpp<1> >;
      extern template class BinaryChiDerivative<2, Cpp<2> >;
      extern template class BinaryChiDerivative<3, Cpp<3> >;
   }
}
#endif
