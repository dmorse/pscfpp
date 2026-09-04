/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/CPT.h>
#include <rp/fts/analyzer/AverageListAnalyzer.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageListAnalyzer<1,CPT>;
      template class AverageListAnalyzer<2,CPT>;
      template class AverageListAnalyzer<3,CPT>;
   }
}
