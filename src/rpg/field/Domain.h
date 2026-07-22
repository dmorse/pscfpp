#ifndef RPG_DOMAIN_H
#define RPG_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Domain.h>     // base class template
#include <pscf/backends/CUT.h>    // base class template parameter

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Domain<1,CUT>;
      extern template class Domain<2,CUT>;
      extern template class Domain<3,CUT>;
   }
}
#endif
