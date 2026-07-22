#ifndef RPC_BD_STEP_FACTORY_H
#define RPC_BD_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStepFactory.h>
#include <pscf/backends/CPT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdStepFactory<1,CPT>;
      extern template class BdStepFactory<2,CPT>;
      extern template class BdStepFactory<3,CPT>;
   }
}
#endif
