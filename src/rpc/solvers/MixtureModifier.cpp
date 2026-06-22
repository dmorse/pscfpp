/*
* PSCF - MixtureModifier Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureModifier.h"
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/solvers/Block.h>

#include <rp/solvers/MixtureModifier.tpp>

namespace Pscf {
   namespace Rp { 
      template class MixtureModifier< Rp::Mixture<1, Rpc::Types<1> > >;
      template class MixtureModifier< Rp::Mixture<2, Rpc::Types<2> > >;
      template class MixtureModifier< Rp::Mixture<3, Rpc::Types<3> > >;
   }
   namespace Rpc { 
      template class MixtureModifier<1>;
      template class MixtureModifier<2>;
      template class MixtureModifier<3>;
   }
}
