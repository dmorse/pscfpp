#ifndef RPG_SCFT_THERMO_H
#define RPG_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/ScftThermo.h>         // class template
#include <pscf/backends/CUT.h>           // template argument
#include <rp/system/SystemConstRef.h>  // base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::ScftThermo<1,CUT>;
      extern template class Rp::ScftThermo<2,CUT>;
      extern template class Rp::ScftThermo<3,CUT>;
   }
} 
#endif
