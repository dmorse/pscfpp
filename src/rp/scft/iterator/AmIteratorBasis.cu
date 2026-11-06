/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/AmIteratorBasis.tpp>
#include <pscf/backend/cuda/CUT.h>

// Explicit instantiation definitions
namespace Pscf {

   // Base class 
   template class 
   AmIteratorTmpl< Rp::Iterator<1,CUT>, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<2,CUT>, DArray<double> >;
   template class 
   AmIteratorTmpl< Rp::Iterator<3,CUT>, DArray<double> >;

   namespace Rp {
      template class AmIteratorBasis<1,CUT>;
      template class AmIteratorBasis<2,CUT>;
      template class AmIteratorBasis<3,CUT>;
   }

}
