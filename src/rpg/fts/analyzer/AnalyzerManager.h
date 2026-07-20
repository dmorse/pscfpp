#ifndef RPG_ANALYZER_MANAGER_H
#define RPG_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <pscf/cuda/CudaTp.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rpg/fts/analyzer/Analyzer.h>       // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, CudaTp<1> >;
      extern template class AnalyzerManager<2, CudaTp<2> >;
      extern template class AnalyzerManager<3, CudaTp<3> >;
   }
}
#endif
