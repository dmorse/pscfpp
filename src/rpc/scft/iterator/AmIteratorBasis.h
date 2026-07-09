#ifndef RPC_AM_ITERATOR_BASIS_H
#define RPC_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.h>  // base class template
#include <rpc/system/Types.h>                  // base class argument 
#include <rpc/scft/iterator/Iterator.h>        // base class member

// Explicit instantiation declarations
namespace Pscf {

   extern template 
   class AmIteratorTmpl<Rp::Iterator<1, Rpc::Types<1> >, DArray<double> >;
   extern template class 
   AmIteratorTmpl<Rp::Iterator<2, Rpc::Types<2> >, DArray<double> >;
   extern template 
   class AmIteratorTmpl<Rp::Iterator<3, Rpc::Types<3> >, DArray<double> >;

   namespace Rp {
      extern template class AmIteratorBasis<1, Rpc::Types<1> >;
      extern template class AmIteratorBasis<2, Rpc::Types<2> >;
      extern template class AmIteratorBasis<3, Rpc::Types<3> >;
   }

}
#endif
