/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

//#include <rpc/fts/brownian/BdStep.h>
//#include <rp/fts/simulator/Simulator.h>
//#include <rp/fts/brownian/BdStepFactory.h>

#include <rpc/fts/analyzer/AnalyzerManager.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>
#include <rpc/fts/trajectory/TrajectoryReader.h>
#include <rpc/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/ramp/Ramp.h>

#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/brownian/BdSimulator.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdSimulator<1,CPT>;
      template class BdSimulator<2,CPT>;
      template class BdSimulator<3,CPT>;
   }
}
