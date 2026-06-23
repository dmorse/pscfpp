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

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Default constructor.
   */
   template <int D>
   SweepParameter<D>::SweepParameter()
    : Rp::SweepParameter< D, Types<D> >()
   {}

   /*
   * Constructor, creates association with system.
   */
   template <int D>
   SweepParameter<D>::SweepParameter(Rp::System<D, Rpc::Types<D> >& system)
    : Rp::SweepParameter< D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class SweepParameter<1, Rpc::Types<1> >;
      template class SweepParameter<2, Rpc::Types<2> >;
      template class SweepParameter<3, Rpc::Types<3> >;
      
   }
   namespace Rpc {
      template class SweepParameter<1>;
      template class SweepParameter<2>;
      template class SweepParameter<3>;
   }
}
