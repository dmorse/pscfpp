/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf {

   template class 
   BlockTmpl< Rp::Propagator<1, CudaTp<1> >, Prdc::RField<1, CudaTp<1> > >;
   template class 
   BlockTmpl< Rp::Propagator<2, CudaTp<2> >, Prdc::RField<2, CudaTp<2> > >;
   template class 
   BlockTmpl< Rp::Propagator<3, CudaTp<3> >, Prdc::RField<3, CudaTp<3> > >;

   namespace Rp {
      template class Block<1, CudaTp<1> >;
      template class Block<2, CudaTp<2> >;
      template class Block<3, CudaTp<3> >;
   }

}
