/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <rp/scft/iterator/AmIteratorBasis.tpp>     // base class implementation

#if 0
// Constructor definition
namespace Pscf {
   namespace Rpg {

      // Constructor
      template <int D>
      Rp::AmIteratorBasis<D, Rpg::Types<D> >::AmIteratorBasis(Rp::System<D, Rpg::Types<D> >& system)
       : Rp::AmIteratorBasis<D, Types<D> >(system)
      {}

   }
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   template class AmIteratorTmpl< Rp::Iterator<1, Rpg::Types<1> >, DArray<double> >;
   template class AmIteratorTmpl< Rp::Iterator<2, Rpg::Types<2> >, DArray<double> >;
   template class AmIteratorTmpl< Rp::Iterator<3, Rpg::Types<3> >, DArray<double> >;
   namespace Rp {
      template class AmIteratorBasis<1, Rpg::Types<1> >;
      template class AmIteratorBasis<2, Rpg::Types<2> >;
      template class AmIteratorBasis<3, Rpg::Types<3> >;
   }
}
