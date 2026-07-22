#ifndef RPG_ITERATOR_H
#define RPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>  // base class template
#include <pscf/backends/CUT.h>           // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Iterator<1,CUT>;
      extern template class Iterator<2,CUT>;
      extern template class Iterator<3,CUT>;
   }
} 
#endif
