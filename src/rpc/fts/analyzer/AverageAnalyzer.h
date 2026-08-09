#ifndef RPC_AVERAGE_ANALYZER_H
#define RPC_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h>  // base class template
#include <pscf/backends/CPT.h>                 // base class argument
#include <rp/fts/analyzer/Analyzer.h>        // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageAnalyzer<1,CPT>;
      extern template class AverageAnalyzer<2,CPT>;
      extern template class AverageAnalyzer<3,CPT>;
   }
}
#endif
