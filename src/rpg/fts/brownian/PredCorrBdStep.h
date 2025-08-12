#ifndef RPG_PRED_CORR_BD_STEP_H
#define RPG_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"

#include <prdc/cuda/RField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Predictor-corrector Brownian dynamics stepper.
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class PredCorrBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~PredCorrBdStep();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Setup before simulation.
      */
      virtual void setup();

      /**
      * Take a single Brownian dynamics step.
      * 
      * \return true if converged, false if failed to converge.
      */
      virtual bool step();

   protected:

      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::cudaRandom;
      using ParamComposite::read;

   private:

      // Predictor value of fields (monomer fields)
      DArray< RField<D> > wp_;

      // Correctd (new) values of fields (monomer fields)
      DArray< RField<D> > wf_;

      // Initial deterministic forces (components)
      DArray< RField<D> > dci_;

      // Random displacement components (components)
      DArray< RField<D> > eta_;

      // Change in one component of wc 
      RField<D> dwc_;

      // Change in pressure field component 
      RField<D> dwp_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

   #ifndef RPG_PRED_CORR_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class PredCorrBdStep<1>;
   extern template class PredCorrBdStep<2>;
   extern template class PredCorrBdStep<3>;
   #endif

}
}
#endif
