#ifndef RPC_HAMILTONIAN_ANALYZER_H
#define RPC_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/HamiltonianAnalyzer.h>  // base class template
#include <pscf/cpu/Cpp.h>                     // template argument
#include "AverageListAnalyzer.h"                  // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class HamiltonianAnalyzer< 1, Cpp<1> >;
      extern template class HamiltonianAnalyzer< 2, Cpp<2> >;
      extern template class HamiltonianAnalyzer< 3, Cpp<3> >;
   }
}
#endif
