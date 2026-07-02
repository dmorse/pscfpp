/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"                    // header
#include <rpg/system/System.h>
#include <rp/fts/analyzer/AverageListAnalyzer.tpp>  // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageListAnalyzer<1, Rpg::Types<1> >;
      template class AverageListAnalyzer<2, Rpg::Types<2> >;
      template class AverageListAnalyzer<3, Rpg::Types<3> >;
   }
}
