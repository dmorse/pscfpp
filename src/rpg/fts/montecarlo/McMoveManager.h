#ifndef RPG_MC_MOVE_MANAGER_H
#define RPG_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMoveManager.h> // base class template
#include <rpg/system/Types.h>                // base class template argument
#include <util/param/Manager.h>              // indirect base class
#include <util/containers/DArray.h>          // member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMoveManager<1, Rpg::Types<1> >;
      extern template class McMoveManager<2, Rpg::Types<2> >;
      extern template class McMoveManager<3, Rpg::Types<3> >;
   }
}
#endif
