#ifndef RPG_ANALYZER_H
#define RPG_ANALYZER_H

#include <rp/fts/analyzer/Analyzer.h>   // class template
#include <rpg/system/Types.h>           // class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Analyzer<1, Rpg::Types<1> >;
      extern template class Analyzer<2, Rpg::Types<2> >;
      extern template class Analyzer<3, Rpg::Types<3> >;
   }
}
#endif
