#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>         // class template
#include <rpg/system/Types.h>           // template argument
#include <rpg/system/SystemConstRef.h>  // base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ScftThermo<1, Rpg::Types<1> >;
      extern template class Rp::ScftThermo<2, Rpg::Types<2> >;
      extern template class Rp::ScftThermo<3, Rpg::Types<3> >;
   }
} 
#endif
