/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/Analyzer.tpp>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class Analyzer<1, Rpg::Types<1> >;
      template class Analyzer<2, Rpg::Types<2> >;
      template class Analyzer<3, Rpg::Types<3> >;
   }
}
