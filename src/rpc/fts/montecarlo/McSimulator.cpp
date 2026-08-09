/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveFactory.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/fts/analyzer/Analyzer.h>
#include <rp/fts/analyzer/AnalyzerFactory.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/trajectory/TrajectoryReaderFactory.h>
#include <rp/fts/trajectory/TrajectoryReader.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/simulator/SimState.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/fts/ramp/Ramp.h>
#include <rp/system/System.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/montecarlo/McSimulator.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McSimulator<1,CPT>;
      template class McSimulator<2,CPT>;
      template class McSimulator<3,CPT>;
   }
}
