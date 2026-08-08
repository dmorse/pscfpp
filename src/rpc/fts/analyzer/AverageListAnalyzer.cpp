/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/analyzer/AverageListAnalyzer.h>
#include <rp/system/System.h>

#include <rp/fts/analyzer/AverageListAnalyzer.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageListAnalyzer<1,CPT>;
      template class AverageListAnalyzer<2,CPT>;
      template class AverageListAnalyzer<3,CPT>;
   }
}
