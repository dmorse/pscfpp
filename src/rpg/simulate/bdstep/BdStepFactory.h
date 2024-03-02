#ifndef RPG_BD_STEP_FACTORY_H
#define RPG_BD_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/simulate/bdstep/BdStep.h>

#include <string>

namespace Pscf {
namespace Rpg {

   template <int D> class BdSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of BdStep.
   *
   * \ingroup Rpg_Simulate_BdStep_Module
   */
   template <int D>
   class BdStepFactory : public Factory< BdStep<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent BdSimulator<D> object
      */
      BdStepFactory(BdSimulator<D>& simulator);

      /**
      * Method to create any BdStep supplied with PSCF.
      *
      * \param className name of the BdStep subclass
      * \return BdStep* pointer to new instance of className
      */
      BdStep<D>* factory(const std::string &className) const;

      using Factory< BdStep<D> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      BdSimulator<D>* simulatorPtr_;

   };

   #ifndef RPG_BD_STEP_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class BdStepFactory<1>;
   extern template class BdStepFactory<2>;
   extern template class BdStepFactory<3>;
   #endif

}
}
#endif
