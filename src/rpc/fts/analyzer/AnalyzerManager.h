#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <pscf/cpu/CppTp.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rpc/fts/analyzer/Analyzer.h>       // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, CppTp<1> >;
      extern template class AnalyzerManager<2, CppTp<2> >;
      extern template class AnalyzerManager<3, CppTp<3> >;
   }
}
#endif
