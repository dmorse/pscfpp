#ifndef RPC_MC_MOVE_MANAGER_H
#define RPC_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMoveManager.h> // base class template
#include <pscf/backends/CPT.h>                // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMoveManager<1,CPT>;
      extern template class McMoveManager<2,CPT>;
      extern template class McMoveManager<3,CPT>;
   }
}
#endif
