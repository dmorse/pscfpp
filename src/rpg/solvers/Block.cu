/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf {

   template class 
   BlockTmpl< Rp::Propagator<1,CUT>, Prdc::RField<1,CUT> >;
   template class 
   BlockTmpl< Rp::Propagator<2,CUT>, Prdc::RField<2,CUT> >;
   template class 
   BlockTmpl< Rp::Propagator<3,CUT>, Prdc::RField<3,CUT> >;

   namespace Rp {
      template class Block<1,CUT>;
      template class Block<2,CUT>;
      template class Block<3,CUT>;
   }

}
