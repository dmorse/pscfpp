#ifndef RPG_HAMILTONIAN_ANALYZER_H
#define RPG_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"                  // indirect base class
#include <rp/fts/analyzer/HamiltonianAnalyzer.h>  // base class template
#include <pscf/backends/CUT.h>                     // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class HamiltonianAnalyzer<1,CUT>;
      extern template class HamiltonianAnalyzer<2,CUT>;
      extern template class HamiltonianAnalyzer<3,CUT>;
   }
}
#endif
