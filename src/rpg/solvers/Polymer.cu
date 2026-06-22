/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rpg/solvers/Block.h>
#include <rpg/solvers/Propagator.h>
#include <prdc/cuda/RField.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>

#include <pscf/solvers/PolymerTmpl.tpp>
#include <rp/solvers/Polymer.tpp>

namespace Pscf {
   template class PolymerTmpl< Rpg::Block<1>, Rp::Propagator<1, Rpg::Types<1> > >;
   template class PolymerTmpl< Rpg::Block<2>, Rp::Propagator<2, Rpg::Types<2> > >;
   template class PolymerTmpl< Rpg::Block<3>, Rp::Propagator<3, Rpg::Types<3> > >;
   namespace Rp {
      template class Polymer< 1, Rpg::Types<1> >;
      template class Polymer< 2, Rpg::Types<2> >;
      template class Polymer< 3, Rpg::Types<3> >;
   }
}
