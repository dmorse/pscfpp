#ifndef RPC_ANALYZER_H
#define RPC_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>  // class template
#include <pscf/cpu/CppTp.h>          // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, CppTp<1> >;
      extern template class Analyzer<2, CppTp<2> >;
      extern template class Analyzer<3, CppTp<3> >;
   }
}
#endif
