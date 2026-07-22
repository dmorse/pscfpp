#ifndef RPG_SOLVENT_H
#define RPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.h>    // base class template
#include <pscf/backends/CUT.h>      // base class template parameter
#include <prdc/field/cuda/RField.h>      // member of base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Solvent<1,CUT>;
      extern template class Solvent<2,CUT>;
      extern template class Solvent<3,CUT>;
   }
}
#endif
