/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include <rpg/solvers/Block.h>
#include <rpg/solvers/Propagator.h>
#include <prdc/field/cuda/RField.h>
#include <pscf/solvers/PolymerTmpl.tpp>
#include <pscf/chem/PolymerModel.h>

#include <pscf/solvers/PolymerTmpl.tpp>
#include <rp/solvers/Polymer.tpp>

namespace Pscf {
   template 
   class PolymerTmpl< Rp::Block<1,CUT>, Rp::Propagator<1,CUT> >;
   template 
   class PolymerTmpl< Rp::Block<2,CUT>, Rp::Propagator<2,CUT> >;
   template 
   class PolymerTmpl< Rp::Block<3,CUT>, Rp::Propagator<3,CUT> >;
   namespace Rp {
      template class Polymer<1,CUT>;
      template class Polymer<2,CUT>;
      template class Polymer<3,CUT>;
   }
}
