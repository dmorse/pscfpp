#ifndef RPC_ANALYZER_H
#define RPC_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>  // class template
#include <pscf/backends/CPT.h>          // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1,CPT>;
      extern template class Analyzer<2,CPT>;
      extern template class Analyzer<3,CPT>;
   }
}
#endif
