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

#include <rp/scft/iterator/AmIteratorBasis.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {

   // Base class instantiation definitions
   template class 
   AmIteratorTmpl< Rp::Iterator<1,CPT>, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2,CPT>, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3,CPT>, DArray<double> >;

   namespace Rp {
      template class AmIteratorBasis<1,CPT>;
      template class AmIteratorBasis<2,CPT>;
      template class AmIteratorBasis<3,CPT>;
   }

}
