/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RampParameter.h"
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/solvers/Block.h>
#include <rpc/field/Domain.h>

#include <rp/fts/ramp/RampParameter.tpp>    // base class implementation

namespace Pscf {
namespace Rpc {

   // Default constructor.
   template <int D>
   RampParameter<D>::RampParameter()
    : Rp::RampParameter<D, Types<D> >()
   {}

   // Constructor - creates association with a Simulator.
   template <int D>
   RampParameter<D>::RampParameter(Rp::Simulator<D, Rpc::Types<D> >& simulator)
    : Rp::RampParameter<D, Types<D> >(simulator)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class RampParameter<1, Rpc::Types<1> >;
      template class RampParameter<2, Rpc::Types<2> >;
      template class RampParameter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class RampParameter<1>;
      template class RampParameter<2>;
      template class RampParameter<3>;
   }
}
