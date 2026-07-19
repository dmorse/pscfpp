#ifndef RPC_MC_SIMULATOR_H
#define RPC_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McSimulator.h>      // base class template
#include <pscf/cpu/Cpp.h>                   // template argument
#include <rpc/fts/simulator/Simulator.h>        // indirect base class
#include <rpc/fts/montecarlo/McMoveManager.h>   // member of base class
#include <rpc/fts/analyzer/AnalyzerManager.h>   // member of base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McSimulator<1, Cpp<1> >;
      extern template class McSimulator<2, Cpp<2> >;
      extern template class McSimulator<3, Cpp<3> >;
   }
}
#endif
