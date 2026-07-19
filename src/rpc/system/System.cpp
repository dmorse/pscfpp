/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <rpc/environment/EnvironmentFactory.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/simulator/SimulatorFactory.h>
#include <rpc/fts/compressor/CompressorFactory.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/scft/sweep/Sweep.h>
#include <rpc/scft/sweep/SweepFactory.h>
#include <rpc/scft/iterator/Iterator.h>
#include <rpc/scft/iterator/IteratorFactory.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/Mask.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <prdc/field/cpu/WaveList.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldDft.h>

#include <rp/system/System.tpp>  // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class System< 1, Cpp<1> >;
      template class System< 2, Cpp<2> >;
      template class System< 3, Cpp<3> >;
   }
}
