/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf { 

   template class BlockTmpl< Rpc::Propagator<1>, Prdc::Cpu::RField<1> >;
   template class BlockTmpl< Rpc::Propagator<2>, Prdc::Cpu::RField<2> >;
   template class BlockTmpl< Rpc::Propagator<3>, Prdc::Cpu::RField<3> >;

   namespace Rpc {
      template class Block<1>;
      template class Block<2>;
      template class Block<3>;
   }
}
