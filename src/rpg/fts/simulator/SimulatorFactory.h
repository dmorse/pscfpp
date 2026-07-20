#ifndef RPG_SIMULATOR_FACTORY_H
#define RPG_SIMULATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/SimulatorFactory.h>
#include <rpg/fts/simulator/Simulator.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SimulatorFactory<1, CudaTp<1> >;
      extern template class SimulatorFactory<2, CudaTp<2> >;
      extern template class SimulatorFactory<3, CudaTp<3> >;
   }
}
#endif
