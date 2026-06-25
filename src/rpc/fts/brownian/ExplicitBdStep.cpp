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

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   ExplicitBdStep<D>::ExplicitBdStep(Rp::BdSimulator<D, Rpc::Types<D> >& simulator)
    : Rp::ExplicitBdStep<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Rp::ExplicitBdStep<1, Rpc::Types<1> >;
      template class Rp::ExplicitBdStep<2, Rpc::Types<2> >;
      template class Rp::ExplicitBdStep<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ExplicitBdStep<1>;
      template class ExplicitBdStep<2>;
      template class ExplicitBdStep<3>;
   }
}
