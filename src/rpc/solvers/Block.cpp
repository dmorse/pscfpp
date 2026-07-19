/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf { 

   template class BlockTmpl< Rp::Propagator<1, Cpp<1> >, Prdc::Cpu::RField<1> >;
   template class BlockTmpl< Rp::Propagator<2, Cpp<2> >, Prdc::Cpu::RField<2> >;
   template class BlockTmpl< Rp::Propagator<3, Cpp<3> >, Prdc::Cpu::RField<3> >;

   namespace Rp {
      template class Block<1, Cpp<1> >;
      template class Block<2, Cpp<2> >;
      template class Block<3, Cpp<3> >;
   }
}
