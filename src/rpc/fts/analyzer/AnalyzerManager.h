#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <rpc/system/Types.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rpc/fts/analyzer/Analyzer.h>       // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, Rpc::Types<1> >;
      extern template class AnalyzerManager<2, Rpc::Types<2> >;
      extern template class AnalyzerManager<3, Rpc::Types<3> >;
   }
}
#endif
