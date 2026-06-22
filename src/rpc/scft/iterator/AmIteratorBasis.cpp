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
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <rpc/field/Mask.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/iterator/AmIteratorBasis.tpp>     // base class implementation

namespace Pscf {
   namespace Rpc {
  
      // Constructor
      template <int D>
      AmIteratorBasis<D>::AmIteratorBasis(System<D>& system)
       : Rp::AmIteratorBasis<D, Types<D> >(system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rp::Iterator<1, Rpc::Types<1> >, DArray<double> >;
   template class AmIteratorTmpl< Rp::Iterator<2, Rpc::Types<2> >, DArray<double> >;
   template class AmIteratorTmpl< Rp::Iterator<3, Rpc::Types<3> >, DArray<double> >;
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
