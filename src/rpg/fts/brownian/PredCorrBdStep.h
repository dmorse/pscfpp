#ifndef RPG_PRED_CORR_BD_STEP_H
#define RPG_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/PredCorrBdStep.h> // base class template
#include <rpg/system/Types.h>               // base class template argument 
#include <prdc/cuda/RField.h>               // base class member
#include <rpg/fts/brownian/BdStep.h>        // indirect base class

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class BdSimulator;

   /**
   * Predictor-corrector Brownian dynamics time stepper.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::PredCorrBdStep, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::PredCorrBdStep
   * \see \ref rp_PredCorrBdStep_page "Manual Page"
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class PredCorrBdStep : public Rp::PredCorrBdStep<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(BdSimulator<D>& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Rp::PredCorrBdStep<1, Rpg::Types<1> >;
      extern template class Rp::PredCorrBdStep<2, Rpg::Types<2> >;
      extern template class Rp::PredCorrBdStep<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class PredCorrBdStep<1>;
      extern template class PredCorrBdStep<2>;
      extern template class PredCorrBdStep<3>;
   }
}
#endif
