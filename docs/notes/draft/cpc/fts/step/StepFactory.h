#ifndef CPC_STEP_FACTORY_H
#define CPC_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <cpc/fts/step/Step.h>

#include <string>

namespace Pscf {
namespace Cpc {

   template <int D> class Simulator;

   using namespace Util;

   /**
   * Factory for subclasses of Step.
   *
   * \ingroup Cpc_Fts_Step_Module
   */
   template <int D>
   class StepFactory : public Factory< Step<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Simulator<D> object
      */
      StepFactory(Simulator<D>& simulator);

      /**
      * Method to create any Step supplied with PSCF.
      *
      * \param className name of the Step subclass
      * \return Step* pointer to new instance of className
      */
      Step<D>* factory(const std::string &className) const;

      using Factory< Step<D> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Simulator<D>* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class StepFactory<1>;
   extern template class StepFactory<2>;
   extern template class StepFactory<3>;

}
}
#endif
