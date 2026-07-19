/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"                   // class header

#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <pscf/cpu/CpuVecRandom.h>
#include <pscf/cpu/VecOp.h>

#include <rp/fts/brownian/ExplicitBdStep.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ExplicitBdStep<1, Cpp<1> >;
      template class ExplicitBdStep<2, Cpp<2> >;
      template class ExplicitBdStep<3, Cpp<3> >;
   }
}
