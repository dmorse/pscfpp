#ifndef RPG_ANALYZER_MANAGER_H
#define RPG_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // class template
#include <rpg/system/Types.h>                // template argument
#include <util/param/Manager.h>              // base class template
#include <rpg/fts/analyzer/Analyzer.h>       // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, Rpg::Types<1> >;
      extern template class AnalyzerManager<2, Rpg::Types<2> >;
      extern template class AnalyzerManager<3, Rpg::Types<3> >;
   }
}
#endif
