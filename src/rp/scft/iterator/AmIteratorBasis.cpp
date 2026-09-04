/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>

#include <pscf/backend/cpp/CPT.h>
#include <rp/scft/iterator/AmIteratorBasis.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {

   // Base class 
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
