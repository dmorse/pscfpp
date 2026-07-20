#ifndef RPG_PERTURBATION_DERIVATIVE_H
#define RPG_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rp/fts/analyzer/PerturbationDerivative.h>
#include <pscf/cuda/CudaTp.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationDerivative< 1, CudaTp<1> >;
      extern template class PerturbationDerivative< 2, CudaTp<2> >;
      extern template class PerturbationDerivative< 3, CudaTp<3> >;
   }
}
#endif
