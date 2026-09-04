#ifndef RPC_FORCE_BIAS_MOVE_H
#define RPC_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMoveBase.h>   // base class template
#include <pscf/backend/CPT.h>           

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   // Declare primary template
   template <int D, class T> class ForceBiasMove;

   /**
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler-Maruyama Brownian
   * dynamics step. 
   *
   * \see ForceBiasMove
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove<D,CPT> 
     : public ForceBiasMoveBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(McSimulator<D,CPT>& simulator);

   private:

      /**
      * Compute force bias field.
      */
      void computeForceBias(RField<D,CPT>& result,
                            RField<D,CPT> const & di,
                            RField<D,CPT> const & df,
                            RField<D,CPT> const & dwc,
                            double mobility);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ForceBiasMoveBase<1,CPT>;
      extern template class ForceBiasMoveBase<2,CPT>;
      extern template class ForceBiasMoveBase<3,CPT>;
      extern template class ForceBiasMove<1,CPT>;
      extern template class ForceBiasMove<2,CPT>;
      extern template class ForceBiasMove<3,CPT>;
   }
}
#endif
