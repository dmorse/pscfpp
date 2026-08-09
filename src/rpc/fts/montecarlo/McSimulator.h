#ifndef RPC_MC_SIMULATOR_H
#define RPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McSimulator.h>      // base class template
#include <pscf/backends/CPT.h>                   // template argument
#include <rp/fts/simulator/Simulator.h>        // indirect base class
#include <rp/fts/montecarlo/McMoveManager.h>   // member of base class
#include <rp/fts/analyzer/AnalyzerManager.h>   // member of base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McSimulator<1,CPT>;
      extern template class McSimulator<2,CPT>;
      extern template class McSimulator<3,CPT>;
   }
}
#endif
