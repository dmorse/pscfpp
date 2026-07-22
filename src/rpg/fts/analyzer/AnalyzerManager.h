#ifndef RPG_ANALYZER_MANAGER_H
#define RPG_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <pscf/backends/CUT.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rpg/fts/analyzer/Analyzer.h>       // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1,CUT>;
      extern template class AnalyzerManager<2,CUT>;
      extern template class AnalyzerManager<3,CUT>;
   }
}
#endif
