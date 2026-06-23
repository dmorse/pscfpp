/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"

#include <rpc/fts/montecarlo/McMove.h>
#include <rpc/fts/montecarlo/McMoveFactory.h>
#include <rpc/fts/montecarlo/McMoveManager.h>
#include <rpc/fts/analyzer/Analyzer.h>
#include <rpc/fts/analyzer/AnalyzerFactory.h>
#include <rpc/fts/analyzer/AnalyzerManager.h>
#include <rpc/fts/trajectory/TrajectoryReaderFactory.h>
#include <rpc/fts/trajectory/TrajectoryReader.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/system/System.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <prdc/cpu/RField.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/montecarlo/McSimulator.tpp>

namespace Pscf {
   namespace Rpc {
   
      /*
      * Constructor.
      */
      template <int D>
      McSimulator<D>::McSimulator(Rp::System<D, Rpc::Types<D> >& system)
       : Rp::McSimulator<D, Types<D> >(system, *this)
      {}
   
   }
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class McSimulator< 1, Rpc::Types<1> >;
      template class McSimulator< 2, Rpc::Types<2> >;
      template class McSimulator< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class McSimulator<1>;
      template class McSimulator<2>;
      template class McSimulator<3>;
   }
}
