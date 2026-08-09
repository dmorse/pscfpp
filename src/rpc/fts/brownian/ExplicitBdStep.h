#ifndef RPC_EXPLICIT_BD_STEP_H
#define RPC_EXPLICIT_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/ExplicitBdStep.h>  // base class template
#include <pscf/backends/CPT.h>                // base template argument 
#include <rp/fts/brownian/BdStep.h>         // indirect base class
#include <prdc/field/cpu/RField.h>           // base class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ExplicitBdStep<1,CPT>;
      extern template class ExplicitBdStep<2,CPT>;
      extern template class ExplicitBdStep<3,CPT>;
   }
}
#endif
