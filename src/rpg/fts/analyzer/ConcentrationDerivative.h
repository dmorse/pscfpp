#ifndef RPG_CONCENTRATION_DERIVATIVE_H
#define RPG_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationDerivative.h> // base class template
#include <rpg/system/Types.h>                        // base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>        // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationDerivative<1, Rpg::Types<1> >;
      extern template class ConcentrationDerivative<2, Rpg::Types<2> >;
      extern template class ConcentrationDerivative<3, Rpg::Types<3> >;
   }
}
#endif
