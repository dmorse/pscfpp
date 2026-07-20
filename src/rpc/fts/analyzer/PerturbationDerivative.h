#ifndef RPC_PERTURBATION_DERIVATIVE_H
#define RPC_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/PerturbationDerivative.h> // base class template
#include <pscf/cpu/CppTp.h>                       // base class argument
#include <rp/fts/analyzer/AverageAnalyzer.h>        // indirect base

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationDerivative< 1, CppTp<1> >;
      extern template class PerturbationDerivative< 2, CppTp<2> >;
      extern template class PerturbationDerivative< 3, CppTp<3> >;
   }
}
#endif
