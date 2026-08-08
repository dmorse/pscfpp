/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Subclasses of Rpc::Iterator 
#include <rp/scft/iterator/AmIteratorBasis.h>
#include <rp/scft/iterator/AmIteratorGrid.h>

#include <pscf/backends/CPT.h> // template implementation
#include <rp/scft/iterator/IteratorFactory.tpp> // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IteratorFactory<1,CPT>;
      template class IteratorFactory<2,CPT>;
      template class IteratorFactory<3,CPT>;
   }
}
