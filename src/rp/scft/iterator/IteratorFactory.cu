/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

// Subclasses of Rp::Iterator 
#include <rp/scft/iterator/AmIteratorBasis.h>
#include <rp/scft/iterator/AmIteratorGrid.h>

#include <pscf/backend/cuda/CUT.h>
#include <rp/scft/iterator/IteratorFactory.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IteratorFactory<1,CUT>;
      template class IteratorFactory<2,CUT>;
      template class IteratorFactory<3,CUT>;
   }
}
