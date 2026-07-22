#ifndef RPG_ANALYZER_H
#define RPG_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>   // class template
#include <pscf/backends/CUT.h>           // class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1,CUT>;
      extern template class Analyzer<2,CUT>;
      extern template class Analyzer<3,CUT>;
   }
}
#endif
