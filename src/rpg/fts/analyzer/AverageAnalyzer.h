#ifndef RPG_AVERAGE_ANALYZER_H
#define RPG_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AverageAnalyzer.h> // class template
#include <rpg/system/Types.h>                // template argument
#include <rpg/fts/analyzer/Analyzer.h>       // base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageAnalyzer<1, Rpg::Types<1> >;
      extern template class AverageAnalyzer<2, Rpg::Types<2> >;
      extern template class AverageAnalyzer<3, Rpg::Types<3> >;
   }
}
#endif
