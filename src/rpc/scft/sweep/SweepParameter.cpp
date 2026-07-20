/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepParameter.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/solvers/Block.h>
#include <rpc/field/Domain.h>

#include <rp/scft/sweep/SweepParameter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1, CppTp<1> >;
      template class SweepParameter<2, CppTp<2> >;
      template class SweepParameter<3, CppTp<3> >;
      
   }
}
