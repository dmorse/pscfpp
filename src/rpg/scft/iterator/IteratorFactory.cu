/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IteratorFactory.h"  

// Subclasses of Iterator 
#include <rpg/scft/iterator/AmIteratorBasis.h>
#include <rpg/scft/iterator/AmIteratorGrid.h>

#include <rp/scft/iterator/IteratorFactory.tpp>  // template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IteratorFactory<1, Rpg::Types<1> >;
      template class IteratorFactory<2, Rpg::Types<2> >;
      template class IteratorFactory<3, Rpg::Types<3> >;
   }
}
