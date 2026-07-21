/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"

#include <rpg/fts/simulator/SimState.h>
#include <rpg/fts/compressor/CompressorFactory.h>
#include <rpg/fts/perturbation/Perturbation.h>
#include <rpg/fts/perturbation/PerturbationFactory.h>
#include <rpg/fts/ramp/Ramp.h>
#include <rpg/fts/ramp/RampFactory.h>

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaVecRandom.h>

#include <rp/fts/simulator/Simulator.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Simulator<1, CudaTp<1> >;
      template class Simulator<2, CudaTp<2> >;
      template class Simulator<3, CudaTp<3> >;
   }
}
