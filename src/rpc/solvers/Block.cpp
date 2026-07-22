/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf { 

   template class BlockTmpl< Rp::Propagator<1,CPT>, Prdc::RField<1,CPT> >;
   template class BlockTmpl< Rp::Propagator<2,CPT>, Prdc::RField<2,CPT> >;
   template class BlockTmpl< Rp::Propagator<3,CPT>, Prdc::RField<3,CPT> >;

   namespace Rp {
      template class Block<1,CPT>;
      template class Block<2,CPT>;
      template class Block<3,CPT>;
   }
}
