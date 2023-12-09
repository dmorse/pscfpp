#ifndef PSPC_EULER_BD_STEP_H
#define PSPC_EULER_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"

#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Pspc_Simulate_BdStep_Module
   */
   template <int D>
   class EulerBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param bdSimulator  parent BdSimulator object
      */
      EulerBdStep(BdSimulator<D>& bdSimulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~EulerBdStep();

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
      */
      virtual void step();

   protected:

      using BdStep<D>::system;
      using BdStep<D>::simulator;
      using BdStep<D>::random;
      using ParamComposite::read;

   private:

      // Local copy of w fields
      DArray< RField<D> > w_;

      // Change in one component of wc
      RField<D> dwc_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

   };

   #ifndef PSPC_EULER_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class EulerBdStep<1>;
   extern template class EulerBdStep<2>;
   extern template class EulerBdStep<3>;
   #endif

}
}
#endif
