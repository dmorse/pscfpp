#ifndef RPC_AVERAGE_LIST_ANALYZER_H
#define RPC_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageListAnalyzer.h> // class template
#include <pscf/cpu/Cpp.h>                    // template argument
#include <rpc/fts/analyzer/Analyzer.h>           // base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageListAnalyzer<1, Cpp<1> >;
      extern template class AverageListAnalyzer<2, Cpp<2> >;
      extern template class AverageListAnalyzer<3, Cpp<3> >;
   }
}
#endif
