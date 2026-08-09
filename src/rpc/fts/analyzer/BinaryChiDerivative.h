#ifndef RPC_CHI_DERIVATIVE_H
#define RPC_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryChiDerivative.h>  // base class template
#include <pscf/backends/CPT.h>                     // base template argument
#include <rp/fts/analyzer/AverageAnalyzer.h>     // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryChiDerivative<1,CPT>;
      extern template class BinaryChiDerivative<2,CPT>;
      extern template class BinaryChiDerivative<3,CPT>;
   }
}
#endif
