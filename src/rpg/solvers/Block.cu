/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf {

   template class BlockTmpl< Rp::Propagator<1, Rpg::Types<1> >, Prdc::Cuda::RField<1> >;
   template class BlockTmpl< Rp::Propagator<2, Rpg::Types<2> >, Prdc::Cuda::RField<2> >;
   template class BlockTmpl< Rp::Propagator<3, Rpg::Types<3> >, Prdc::Cuda::RField<3> >;

   namespace Rp {
      template class Block<1, Rpg::Types<1> >;
      template class Block<2, Rpg::Types<2> >;
      template class Block<3, Rpg::Types<3> >;
   }

}
