/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>

#include <rp/fts/analyzer/AverageAnalyzer.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageAnalyzer<1,CUT>;
      template class AverageAnalyzer<2,CUT>;
      template class AverageAnalyzer<3,CUT>;
   }
}
