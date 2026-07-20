/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulatorFactory.h"

// Subclasses of Simulator
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/brownian/BdSimulator.h>

#include <rp/fts/simulator/SimulatorFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SimulatorFactory<1, CudaTp<1> >;
      template class SimulatorFactory<2, CudaTp<2> >;
      template class SimulatorFactory<3, CudaTp<3> >;
   }
}
