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

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Manager for a set of McMove objects.
   *
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class McMoveManager : public Rp::McMoveManager< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      * \param system parent System
      */
      McMoveManager(McSimulator<D>& simulator, Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMoveManager<1, Rpg::Types<1> >;
      extern template class McMoveManager<2, Rpg::Types<2> >;
      extern template class McMoveManager<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class McMoveManager<1>;
      extern template class McMoveManager<2>;
      extern template class McMoveManager<3>;
   }
}
#endif
