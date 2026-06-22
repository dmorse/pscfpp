/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"                    // class header
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <rpc/field/Mask.h>
#include <prdc/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp> // base class implementation

namespace Pscf {
   namespace Rpc {

      // Constructor
      template <int D>
      AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
       : Rp::AmIteratorGrid<D, Types<D> >(system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   template class 
   AmIteratorTmpl< Rp::Iterator<1, Rpc::Types<1> >, DRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2, Rpc::Types<2> >, DRArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3, Rpc::Types<3> >, DRArray<double> >;
   namespace Rp {
      template class AmIteratorGrid<1, Rpc::Types<1> >;
      template class AmIteratorGrid<2, Rpc::Types<2> >;
      template class AmIteratorGrid<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class AmIteratorGrid<1>;
      template class AmIteratorGrid<2>;
      template class AmIteratorGrid<3>;
   }
}
