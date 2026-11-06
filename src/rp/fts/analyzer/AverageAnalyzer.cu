/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <rp/fts/analyzer/AverageAnalyzer.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageAnalyzer<1,CUT>;
      template class AverageAnalyzer<2,CUT>;
      template class AverageAnalyzer<3,CUT>;
   }
}
