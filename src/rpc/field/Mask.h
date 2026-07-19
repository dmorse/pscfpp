#ifndef RPC_MASK_H
#define RPC_MASK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Mask.h>           // base class template
#include <pscf/cpu/Cpp.h>        // base class template argument
#include <prdc/field/cpu/RField.h>   // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Mask< 1, Cpp<1> >;
      extern template class Mask< 2, Cpp<2> >;
      extern template class Mask< 3, Cpp<3> >;
   }
}
#endif
