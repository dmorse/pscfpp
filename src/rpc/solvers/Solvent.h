#ifndef RPC_SOLVENT_H
#define RPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.h>       // base class template
#include <pscf/backends/CPT.h>         // base class argument
#include <prdc/field/cpu/RField.h>    // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Pscf::Prdc;
      extern template class Solvent<1,CPT>;
      extern template class Solvent<2,CPT>;
      extern template class Solvent<3,CPT>;
   }
}
#endif
