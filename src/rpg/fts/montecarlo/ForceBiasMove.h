#ifndef RPG_FORCE_BIAS_MOVE_H
#define RPG_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMoveBase.h> // base class template
#include <rpg/system/Types.h>                    // base class argument
#include <prdc/field/cuda/RField.h>              // base class member
#include <util/containers/DArray.h>              // base base member
#include <rp/fts/montecarlo/McMove.h>            // indirect base class

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * ForceBiasMove attempts a Brownian dynamics move.
   *
   * This class implements a Monte Carlo move in which the unconstrained
   * attempted move is created by an explicit Euler-Maruyama Brownian
   * dynamics step.
   *
   * Because the probability of attempting a move is not equal to that
   * of generating the reverse move, the acceptance criterion used in
   * the move() function must take into account the ratio of generation
   * probabilities.
   *
   * \see ForceBiasMoveBase
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove<D, Rpg::Types<D> >
    : public ForceBiasMoveBase<D, Rpg::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(McSimulator<D, Rpg::Types<D> >& simulator);

   private:

      /**
      * Compute force bias field.
      */
      void computeForceBias(RField<D> & result,
                            RField<D> const & di,
                            RField<D> const & df,
                            RField<D> const & dwc,
                            double mobility);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ForceBiasMoveBase<1, Rpg::Types<1> >;
      extern template class ForceBiasMoveBase<2, Rpg::Types<2> >;
      extern template class ForceBiasMoveBase<3, Rpg::Types<3> >;
      extern template class ForceBiasMove<1, Rpg::Types<1> >;
      extern template class ForceBiasMove<2, Rpg::Types<2> >;
      extern template class ForceBiasMove<3, Rpg::Types<3> >;
   }
}
#endif
