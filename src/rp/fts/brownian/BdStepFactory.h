#ifndef RP_BD_STEP_FACTORY_H
#define RP_BD_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rp/fts/brownian/BdStep.h>

namespace Pscf {
namespace Rp {

   // Forward declaration
   template <int D, class T> class BdSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of BdStep.
   *
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class BdStepFactory : public Factory< BdStep<D,T> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent BdSimulator<D,T> object
      */
      BdStepFactory(BdSimulator<D,T>& simulator);

      /**
      * Method to create any BdStep supplied with PSCF.
      *
      * \param className name of the BdStep subclass
      * \return BdStep* pointer to new instance of className
      */
      BdStep<D,T>* factory(const std::string &className) const;

      using Factory< BdStep<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      BdSimulator<D,T>* simulatorPtr_;

   };

}
}
#endif
