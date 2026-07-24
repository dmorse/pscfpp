#ifndef RPC_SIMULATOR_FACTORY_H
#define RPC_SIMULATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/SimulatorFactory.h>
#include <pscf/backends/CPT.h>
//#include <rpc/fts/simulator/Simulator.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SimulatorFactory<1,CPT>;
      extern template class SimulatorFactory<2,CPT>;
      extern template class SimulatorFactory<3,CPT>;
   }
}
#endif
