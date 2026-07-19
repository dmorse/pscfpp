#ifndef RPG_C_FIELD_CONTAINER_H
#define RPG_C_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>     // class template
#include <rpg/system/Types.h>     // class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template class Rp::CFields<1, Rpg::Types<1> >;
      extern template class Rp::CFields<2, Rpg::Types<2> >;
      extern template class Rp::CFields<3, Rpg::Types<3> >;
   }
}
#endif
