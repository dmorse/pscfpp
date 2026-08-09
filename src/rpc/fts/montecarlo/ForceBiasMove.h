#ifndef RPC_FORCE_BIAS_MOVE_H
#define RPC_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMoveBase.h>   // base class template
#include <pscf/backends/CPT.h>                      // base class argument
#include <rp/fts/montecarlo/McMove.h>             // indirect base class
#include <prdc/field/cpu/RField.h>                 // base class member

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
   * Specializations of this template with D=1, 2, and 3 are derived
   * from specializations of base class template Rp::ForceBiasMove,
   * and inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::ForceBiasMove
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove<D,CPT> 
     : public Rp::ForceBiasMoveBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(Rp::McSimulator<D,CPT>& simulator);

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
