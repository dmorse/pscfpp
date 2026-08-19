#ifndef RP_LM_BD_STEP_H
#define RP_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStep.h>      // base class
#include <prdc/field/RField.h>           // member
#include <util/containers/DArray.h>      // member

#include <pscf/backends/TmplDeclare.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Leimkuhler-Matthews Brownian dynamics stepper.
   *
   * The Leimkuhler-Matthews step algorithm differs from an explicit
   * Euler algorithm in that it uses a random displacement that is
   * given by a sum of random numbers generated at this step and the
   * previous step.
   *
   * References:
   *
   *  - B. Vorselaars, J. Chemical Physics, 158 114117 (2023)
   *    [ https://doi.org/10.1063/5.0131183 ]
   *
   *  - B. Leimkuhler and C. Matthews, Applied Mathematics Research Express,
   *    Issue 1, pages 34-56 (2013) [ https://doi.org/10.1093/amrx/abs010 ]
   *
   *  - B. Leimkuhler and C. Matthews, J. Chemical Physics,
   *    vol. 138, 174102 (2013) [ https://doi.org/10.1063/1.4802990 ]
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : backend identifier class (CPT or CUT)
   *
   * \see \ref rp_LMBdStep_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class LMBdStep : public BdStep<D,T>
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      LMBdStep(BdSimulator<D,T>& simulator);

      ~LMBdStep() = default;

      /**
      * Read body of parameter file block.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream &in) override;

      /**
      * Setup before simulation loop.
      */
      void setup() override;

      /**
      * Take a single Brownian dynamics step.
      *
      * \return true if compressor converged, false otherwise
      */
      bool step() override;

   protected:

      // Alias for base class.
      using BdStepT = BdStep<D,T>;

      // Inherited functions (selected)
      using BdStepT::system;
      using BdStepT::simulator;
      using BdStepT::vecRandom;

   private:

      using RFieldT = RField<D,T>;

      // Private data members

      // New field values
      DArray< RFieldT > w_;

      // Random displacements (A)
      DArray< RFieldT > etaA_;

      // Random displacements (B)
      DArray< RFieldT > etaB_;

      // Change in one field component
      RFieldT dwc_;

      // Pointer to new random displacements
      DArray< RFieldT >* etaNewPtr_;

      // Pointer to old random displacements
      DArray< RFieldT >* etaOldPtr_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

      // Private member functions

      RFieldT& etaNew(int i)
      {  return (*etaNewPtr_)[i]; }

      RFieldT& etaOld(int i)
      {  return (*etaOldPtr_)[i]; }

      /// Generate new values for etaNew
      void generateEtaNew();

      /// Exchange pointer values for etaNew and etaOld.
      void exchangeOldNew();

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(LMBdStep)

}
}
#endif
