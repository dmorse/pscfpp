/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <rpc/fts/simulator/SimState.h>
#include <rpc/fts/compressor/CompressorFactory.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/fts/perturbation/PerturbationFactory.h>
#include <rpc/fts/ramp/Ramp.h>
#include <rpc/fts/ramp/RampFactory.h>

#include <prdc/cpu/RField.h>

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/cpu/CpuVecRandom.h>

#include <rp/fts/simulator/Simulator.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Simulator<1, Rpc::Types<1> >;
      template class Simulator<2, Rpc::Types<2> >;
      template class Simulator<3, Rpc::Types<3> >;
   }
}
