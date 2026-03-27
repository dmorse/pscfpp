#ifndef RP_PRED_CORR_BD_STEP_H
#define RP_PRED_CORR_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/containers/DArray.h>        // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Predictor-corrector Brownian dynamics stepper.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named PredCorrBdStep,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_PredCorrBdStep_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class PredCorrBdStep : public T::BdStep
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      PredCorrBdStep(typename T::BdSimulator& simulator);

      /**
      * Destructor.
      */
      virtual ~PredCorrBdStep();

      /**
      * Read body of parameter file block and initialize.
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

      using BdStepT = typename T::BdStep;

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
