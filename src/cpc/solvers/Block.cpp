/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf { 

   template class BlockTmpl< Cpc::Propagator<1>, Prdc::CField<1, CppTp<1> > >;
   template class BlockTmpl< Cpc::Propagator<2>, Prdc::CField<2, CppTp<2> > >;
   template class BlockTmpl< Cpc::Propagator<3>, Prdc::CField<3, CppTp<3> > >;

   namespace Cpc {
      template class Block<1>;
      template class Block<2>;
      template class Block<3>;
   }

}
