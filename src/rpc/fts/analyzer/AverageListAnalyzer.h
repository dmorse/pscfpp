#ifndef RPC_AVERAGE_LIST_ANALYZER_H
#define RPC_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageListAnalyzer.h> // class template
#include <pscf/backends/CPT.h>                    // template argument
#include <rp/fts/analyzer/Analyzer.h>           // base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageListAnalyzer<1,CPT>;
      extern template class AverageListAnalyzer<2,CPT>;
      extern template class AverageListAnalyzer<3,CPT>;
   }
}
#endif
