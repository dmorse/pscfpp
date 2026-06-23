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
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>

#include <rp/scft/iterator/AmIteratorGrid.tpp>

#if 0
namespace Pscf {
   namespace Rpg {
  
      // Constructor.
      template <int D>
      Rp::AmIteratorGrid<D, Rpg::Types<D> >::AmIteratorGrid(Rp::System<D, Rpg::Types<D> >& system)
       : Rp::AmIteratorGrid<D, Types<D> >(system)
      {}

   }
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   template class 
   AmIteratorTmpl< Rp::Iterator<1, Rpg::Types<1> >, DeviceArray<cudaReal> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2, Rpg::Types<2> >, DeviceArray<cudaReal> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3, Rpg::Types<3> >, DeviceArray<cudaReal> >;
   namespace Rp {
      template class AmIteratorGrid<1, Rpg::Types<1> >;
      template class AmIteratorGrid<2, Rpg::Types<2> >;
      template class AmIteratorGrid<3, Rpg::Types<3> >;
   }
}
