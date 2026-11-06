/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.tpp>
#include <pscf/backend/cpp/CPT.h>

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
