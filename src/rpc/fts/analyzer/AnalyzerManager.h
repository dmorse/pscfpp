#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <pscf/backends/CPT.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rp/fts/analyzer/Analyzer.h>       // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1,CPT>;
      extern template class AnalyzerManager<2,CPT>;
      extern template class AnalyzerManager<3,CPT>;
   }
}
#endif
