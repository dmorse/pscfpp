#ifndef RPC_C_FIELDS_H
#define RPC_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>       // base class template
#include <pscf/cpu/CppTp.h>       // base class parameter

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CFields<1, CppTp<1> >;
      extern template class CFields<2, CppTp<2> >;
      extern template class CFields<3, CppTp<3> >;
   }
}
#endif
