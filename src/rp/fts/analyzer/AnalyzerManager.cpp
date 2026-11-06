/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>
#include <rp/fts/analyzer/AnalyzerManager.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class AnalyzerManager<1,CPT>;
      template class AnalyzerManager<2,CPT>;
      template class AnalyzerManager<3,CPT>;
   }
}
