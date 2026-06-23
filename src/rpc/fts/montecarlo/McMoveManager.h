#ifndef RPC_MC_MOVE_MANAGER_H
#define RPC_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMoveManager.h> // base class template
#include <rpc/system/Types.h>                // base class template argument
#include <util/param/Manager.h>              // indirect base class
#include <util/containers/DArray.h>          // member

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
   namespace Rpc {
      template <int D> class McSimulator;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Manager for a set of McMove objects.
   *
   * \see \ref rp_McMoveManager_page "Manual Page"
   * \ingroup Rpc_Fts_MonteCarlo_Module
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
      McMoveManager(McSimulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMoveManager<1, Rpc::Types<1> >;
      extern template class McMoveManager<2, Rpc::Types<2> >;
      extern template class McMoveManager<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class McMoveManager<1>;
      extern template class McMoveManager<2>;
      extern template class McMoveManager<3>;
   }
}
#endif
