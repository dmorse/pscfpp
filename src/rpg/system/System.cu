#ifndef RPG_SYSTEM_CU
#define RPG_SYSTEM_CU

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"               // class header

#include <rpg/environment/EnvironmentFactory.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/simulator/SimulatorFactory.h>
#include <rpg/scft/ScftThermo.h>
#include <rpg/scft/iterator/Iterator.h>
#include <rpg/scft/iterator/IteratorFactory.h>
#include <rpg/scft/sweep/Sweep.h>
#include <rpg/scft/sweep/SweepFactory.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/MixtureModifier.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>

#include <prdc/field/cuda/RField.h>
#include <prdc/field/cuda/RFieldDft.h>
#include <prdc/field/cuda/WaveList.h>
#include <prdc/environment/Environment.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/ThreadMesh.h>

#include <rp/system/System.tpp>   // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class System< 1, CudaTp<1> >;
      template class System< 2, CudaTp<2> >;
      template class System< 3, CudaTp<3> >;
   }
}
#endif
