/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/iterator/AmIteratorBasis.tpp>     // base class implementation

// Constructor definition
namespace Pscf {
   namespace Rpc {
  
      template <int D>
      AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
       : Rp::AmIteratorBasis<D, Types<D> >(system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rpc::Iterator<1>, DArray<double> >;
   template class AmIteratorTmpl< Rpc::Iterator<2>, DArray<double> >;
   template class AmIteratorTmpl< Rpc::Iterator<3>, DArray<double> >;
   namespace Rp {
      template class AmIteratorBasis<1, Rpc::Types<1> >;
      template class AmIteratorBasis<2, Rpc::Types<2> >;
      template class AmIteratorBasis<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class AmIteratorBasis<1>;
      template class AmIteratorBasis<2>;
      template class AmIteratorBasis<3>;
   }
}
