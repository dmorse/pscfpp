/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.tpp>
#include <pscf/backends/CUT.h>

namespace Pscf {

   // Base class instantiations
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
