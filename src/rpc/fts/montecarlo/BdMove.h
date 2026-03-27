#ifndef RPC_BD_MOVE_H
#define RPC_BD_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/BdMove.h>   // base class template
#include <rpc/system/Types.h>           // base class template argument 
#include <rpc/fts/montecarlo/McMove.h>  // indirect base class
#include <prdc/cpu/RField.h>            // base class member

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class McSimulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Brownian dynamics Monte-Carlo move.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::BdMove < D, Types<D> >, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::BdMove
   * \see \ref rp_BdMove_page "Manual Page"
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class BdMove : public Rp::BdMove<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      BdMove(McSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::BdMove<1, Rpc::Types<1> >;
      extern template class Rp::BdMove<2, Rpc::Types<2> >;
      extern template class Rp::BdMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BdMove<1>;
      extern template class BdMove<2>;
      extern template class BdMove<3>;
   }
}
#endif
