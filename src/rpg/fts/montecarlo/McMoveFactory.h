#ifndef RPG_MC_MOVE_FACTORY_H
#define RPG_MC_MOVE_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMoveFactory.h>
#include <pscf/backends/CUT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMoveFactory<1,CUT>;
      extern template class McMoveFactory<2,CUT>;
      extern template class McMoveFactory<3,CUT>;
   }
}
#endif
