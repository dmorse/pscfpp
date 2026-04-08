#ifndef RPG_FORCE_BIAS_MOVE_H
#define RPG_FORCE_BIAS_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/ForceBiasMove.h> // direct base class template
#include <rpg/system/Types.h>                // direct base argument
#include <prdc/cuda/RField.h>                // direct base member
#include <util/containers/DArray.h>          // direct base member
#include "McMove.h"                          // indirect base class

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class McSimulator;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

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
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from specializations of the base class template Rp::ForceBiasMove, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::ForceBiasMove
   * \see \ref rp_ForceBiasMove_page "Manual Page"
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class ForceBiasMove : public Rp::ForceBiasMove<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator
      */
      ForceBiasMove(McSimulator<D>& simulator);

   private:

      /**
      * Compute force bias field.
      */
      void computeForceBias(RField<D>& result,
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
      extern template class ForceBiasMove<1, Rpg::Types<1> >;
      extern template class ForceBiasMove<2, Rpg::Types<2> >;
      extern template class ForceBiasMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ForceBiasMove<1>;
      extern template class ForceBiasMove<2>;
      extern template class ForceBiasMove<3>;
   }
}
#endif
