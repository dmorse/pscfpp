#ifndef RPG_LM_BD_STEP_H
#define RPG_LM_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "BdStep.h"

#include <prdc/cuda/RField.h>
#include <util/containers/DArray.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Leimkuhler-Matthews Brownian dynamics stepper.
   *
   * As described in:
   * 
   *   B. Vorselaars, J. Chemical Physics, 158 114117 (2023)
   *   [ https://doi.org/10.1063/5.0131183 ]
   *
   *   B. Leimkuhler and C. Matthews, Applied Mathematics Research Express,
   *   Issue 1, pages 34-56 (2013) [ https://doi.org/10.1093/amrx/abs010 ]
   *
   *   B. Leimkuhler and C. Matthews, J. Chemical Physics, 
   *   vol. 138, 174102 (2013) [ https://doi.org/10.1063/1.4802990 ] 
   * 
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class LMBdStep : public BdStep<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      LMBdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~LMBdStep();

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

      // Private data members

      // New field values
      DArray< RField<D> > w_;

      // Random displacements (A)
      DArray< RField<D> > etaA_;

      // Random displacements (B)
      DArray< RField<D> > etaB_;

      // Change in one field component
      RField<D> dwc_;

      // Pointer to new random displacements
      DArray< RField<D> >* etaNewPtr_;

      // Pointer to old random displacements
      DArray< RField<D> >* etaOldPtr_;
      
      // Normal-distributed random fields
      RField<D> gaussianField_;

      // Prefactor of -dc_ in deterministic drift term
      double mobility_;

      // Private member functions

      RField<D>& etaNew(int i) 
      {   return (*etaNewPtr_)[i]; }

      RField<D>& etaOld(int i) 
      {   return (*etaOldPtr_)[i]; }

      /// Generate new values for etaNew
      void generateEtaNew();

      /// Exchange pointer values for etaNew and etaOld.
      void exchangeOldNew();

   };

   #ifndef RPG_LM_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class LMBdStep<1>;
   extern template class LMBdStep<2>;
   extern template class LMBdStep<3>;
   #endif

}
}
#endif
