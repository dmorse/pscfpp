/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/simulator/SimState.h>
#include <rpc/fts/compressor/CompressorFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/perturbation/PerturbationFactory.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/fts/ramp/RampFactory.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/backends/CPT.h>

#include <rp/fts/simulator/Simulator.tpp>  // class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Simulator<1,CPT>;
      template class Simulator<2,CPT>;
      template class Simulator<3,CPT>;
   }
}
