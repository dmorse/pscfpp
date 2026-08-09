/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdSimulator.h"
#include <rp/fts/brownian/BdStep.h>
#include <rp/fts/brownian/BdStepFactory.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/analyzer/AnalyzerFactory.h>
#include <rp/fts/trajectory/TrajectoryReader.h>
#include <rp/fts/trajectory/TrajectoryReaderFactory.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/fts/ramp/Ramp.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

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
