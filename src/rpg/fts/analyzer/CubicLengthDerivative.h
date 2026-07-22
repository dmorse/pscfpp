#ifndef RPG_CUBIC_LENGTH_DERIVATIVE_H
#define RPG_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.h> // direct base template
#include <pscf/backends/CUT.h>                      // direct base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>      // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CubicLengthDerivative<1,CUT>;
      extern template class CubicLengthDerivative<2,CUT>;
      extern template class CubicLengthDerivative<3,CUT>;
   }
}
#endif
