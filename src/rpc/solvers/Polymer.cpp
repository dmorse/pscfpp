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

   template class 
   PolymerTmpl< Rp::Block<1,CPT>, Rp::Propagator<1,CPT> >;
   template class 
   PolymerTmpl< Rp::Block<2,CPT>, Rp::Propagator<2,CPT> >;
   template class 
   PolymerTmpl< Rp::Block<3,CPT>, Rp::Propagator<3,CPT> >;

   namespace Rp {
      template class Polymer<1,CPT>;
      template class Polymer<2,CPT>;
      template class Polymer<3,CPT>;
   }
}
