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

#include <rp/scft/iterator/AmIteratorBasis.tpp>     // implementation

// Explicit instantiation definitions
namespace Pscf {
   template class 
   AmIteratorTmpl< Rp::Iterator<1, CudaTp<1> >, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2, CudaTp<2> >, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3, CudaTp<3> >, DArray<double> >;
   namespace Rp {
      template class AmIteratorBasis<1, CudaTp<1> >;
      template class AmIteratorBasis<2, CudaTp<2> >;
      template class AmIteratorBasis<3, CudaTp<3> >;
   }
}
