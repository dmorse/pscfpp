/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Subclasses of Simulator
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/fts/brownian/BdSimulator.h>

#include <rp/fts/simulator/SimulatorFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SimulatorFactory<1,CPT>;
      template class SimulatorFactory<2,CPT>;
      template class SimulatorFactory<3,CPT>;
   }
}
