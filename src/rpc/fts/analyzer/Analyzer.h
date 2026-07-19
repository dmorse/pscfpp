#ifndef RPC_ANALYZER_H
#define RPC_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>  // class template
#include <pscf/cpu/Cpp.h>          // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, Cpp<1> >;
      extern template class Analyzer<2, Cpp<2> >;
      extern template class Analyzer<3, Cpp<3> >;
   }
}
#endif
