/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf {

   template class BlockTmpl< Rpg::Propagator<1>, Prdc::Cuda::RField<1> >;
   template class BlockTmpl< Rpg::Propagator<2>, Prdc::Cuda::RField<2> >;
   template class BlockTmpl< Rpg::Propagator<3>, Prdc::Cuda::RField<3> >;

   namespace Rpg {
      template class Block<1>;
      template class Block<2>;
      template class Block<3>;
   }

}
