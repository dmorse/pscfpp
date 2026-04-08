#ifndef RPG_BD_MOVE_H
#define RPG_BD_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/BdMove.h>   // base class template
#include <rpg/system/Types.h>           // base class template argument 
#include <rpg/fts/montecarlo/McMove.h>  // indirect base class
#include <prdc/cuda/RField.h>           // base class member

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class McSimulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * Brownian dynamics Monte-Carlo move.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base class template Rp::BdMove, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::BdMove
   * \see \ref rp_BdMove_page "Manual Page"
   * \ingroup Rpg_Fts_MonteCarlo_Module
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
      extern template class Rp::BdMove<1, Rpg::Types<1> >;
      extern template class Rp::BdMove<2, Rpg::Types<2> >;
      extern template class Rp::BdMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BdMove<1>;
      extern template class BdMove<2>;
      extern template class BdMove<3>;
   }
}
#endif
