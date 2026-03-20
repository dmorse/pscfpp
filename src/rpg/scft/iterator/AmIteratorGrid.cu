/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp>

// Constructor definition
namespace Pscf {
   namespace Rpg {
  
      template <int D>
      AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
       : Rp::AmIteratorGrid<D, Types<D> >(system)
      {}

   }
}

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rpg::Iterator<1>, DeviceArray<cudaReal> >;
   template class AmIteratorTmpl< Rpg::Iterator<2>, DeviceArray<cudaReal> >;
   template class AmIteratorTmpl< Rpg::Iterator<3>, DeviceArray<cudaReal> >;
   namespace Rp {
      template class AmIteratorGrid<1, Rpg::Types<1> >;
      template class AmIteratorGrid<2, Rpg::Types<2> >;
      template class AmIteratorGrid<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class AmIteratorGrid<1>;
      template class AmIteratorGrid<2>;
      template class AmIteratorGrid<3>;
   }
}
