/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/Analyzer.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Analyzer<1, Cpp<1> >;
      template class Analyzer<2, Cpp<2> >;
      template class Analyzer<3, Cpp<3> >;
   }
}
