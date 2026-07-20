#ifndef RP_PRED_CORR_BD_STEP_H
#define RP_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>      // base class
#include <util/containers/DArray.h>      // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Predictor-corrector Brownian dynamics stepper.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, CppTp<D> or CudaTp<D>
   *
   * \see \ref rp_PredCorrBdStep_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class PredCorrBdStep : public BdStep<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(BdSimulator<D,T>& simulator);

      ~PredCorrBdStep() = default;

      /**
      * Read body of parameter file block and allocate memory.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true if converged, false if failed to converge
      */
      bool step() override;

   protected:

      // Alias for base class.
      using BdStepT = BdStep<D,T>;

      // Protected inherited member functions
      using BdStepT::system;
      using BdStepT::simulator;

   private:

      using RFieldT = typename T::RField;

      // Predicted values of fields (monomer fields)
      DArray< RFieldT > wp_;

      // Corrected (final) values of fields (monomer fields)
      DArray< RFieldT > wf_;

      // Initial deterministic forces (eigenvector components)
      DArray< RFieldT > dci_;

      // Random displacement components (eigenvector components)
      DArray< RFieldT > eta_;

      // Change in one component of wc
      RFieldT dwc_;

      // Change in pressure field component
      RFieldT dwp_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

}
}
#endif
