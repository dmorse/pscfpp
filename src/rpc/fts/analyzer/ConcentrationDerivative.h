#ifndef RPC_CONCENTRATION_DERIVATIVE_H
#define RPC_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationDerivative.h> // base class template
#include <pscf/cpu/Cpp.h>                        // base class argument
#include <rp/fts/analyzer/AverageAnalyzer.h>         // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationDerivative<1, Cpp<1> >;
      extern template class ConcentrationDerivative<2, Cpp<2> >;
      extern template class ConcentrationDerivative<3, Cpp<3> >;
   }
}
#endif
