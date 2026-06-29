/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rpc/solvers/Block.h>
#include <rpc/solvers/Propagator.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>

#include <pscf/solvers/PolymerTmpl.tpp>
#include <rp/solvers/Polymer.tpp>

namespace Pscf {
   template class PolymerTmpl< Rp::Block<1, Rpc::Types<1> >, Rp::Propagator<1, Rpc::Types<1> > >;
   template class PolymerTmpl< Rp::Block<2, Rpc::Types<2> >, Rp::Propagator<2, Rpc::Types<2> > >;
   template class PolymerTmpl< Rp::Block<3, Rpc::Types<3> >, Rp::Propagator<3, Rpc::Types<3> > >;
   namespace Rp {
      template class Polymer< 1, Rpc::Types<1> >;
      template class Polymer< 2, Rpc::Types<2> >;
      template class Polymer< 3, Rpc::Types<3> >;
   }
}
