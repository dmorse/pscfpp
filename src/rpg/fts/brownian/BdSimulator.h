#ifndef RPG_BD_SIMULATOR_H
#define RPG_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdSimulator.h>       // base class template
#include <pscf/backends/CUT.h>                  // template argument
#include <rpg/fts/simulator/Simulator.h>       // indirect base class
#include <rpg/fts/analyzer/AnalyzerManager.h>  // member of base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdSimulator<1,CUT>;
      extern template class BdSimulator<2,CUT>;
      extern template class BdSimulator<3,CUT>;
   }
}
#endif
