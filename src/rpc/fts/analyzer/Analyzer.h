#ifndef RPC_ANALYZER_H
#define RPC_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>  // class template
#include <rpc/system/Types.h>          // template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, Rpc::Types<1> >;
      extern template class Analyzer<2, Rpc::Types<2> >;
      extern template class Analyzer<3, Rpc::Types<3> >;
   }
}
#endif
