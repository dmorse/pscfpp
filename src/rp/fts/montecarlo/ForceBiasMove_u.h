#ifndef RPG_FORCE_BIAS_MOVE_H
#define RPG_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMoveBase.h> // base class template
#include <pscf/backends/CUT.h>                   // template argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler-Maruyama Brownian
   * dynamics step.
   *
   * \see ForceBiasMoveBase
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove<D,CUT>
    : public ForceBiasMoveBase<D,CUT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(McSimulator<D,CUT>& simulator);

   private:

      /**
      * Compute force bias field.
      */
      void computeForceBias(RField<D,CUT> & result,
                            RField<D,CUT> const & di,
                            RField<D,CUT> const & df,
                            RField<D,CUT> const & dwc,
                            double mobility);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ForceBiasMoveBase<1,CUT>;
      extern template class ForceBiasMoveBase<2,CUT>;
      extern template class ForceBiasMoveBase<3,CUT>;
      extern template class ForceBiasMove<1,CUT>;
      extern template class ForceBiasMove<2,CUT>;
      extern template class ForceBiasMove<3,CUT>;
   }
}
#endif
