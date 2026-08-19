/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/fts/analyzer/HamiltonianAnalyzer.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class HamiltonianAnalyzer<1,CPT>;
      template class HamiltonianAnalyzer<2,CPT>;
      template class HamiltonianAnalyzer<3,CPT>;
   }
}
