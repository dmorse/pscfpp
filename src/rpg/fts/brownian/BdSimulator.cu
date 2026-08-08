/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"

#include <rpg/fts/brownian/BdStep.h>
#include <rpg/fts/brownian/BdStepFactory.h>
#include <rpg/fts/analyzer/AnalyzerManager.h>
#include <rpg/fts/analyzer/AnalyzerFactory.h>
#include <rpg/fts/trajectory/TrajectoryReader.h>
#include <rpg/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/ramp/Ramp.h>
#include <rpg/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/CudaVecRandom.h>

#include <rp/fts/brownian/BdSimulator.tpp>


// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdSimulator<1,CUT>;
      template class BdSimulator<2,CUT>;
      template class BdSimulator<3,CUT>;
   }
}
