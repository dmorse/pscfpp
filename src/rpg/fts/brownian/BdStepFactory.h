#ifndef RPG_BD_STEP_FACTORY_H
#define RPG_BD_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/brownian/BdStep.h>

namespace Pscf {
namespace Rpg {

   template <int D> class BdSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of BdStep.
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class BdStepFactory : public Factory< Rp::BdStep<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Rp::BdSimulator<D, Rpg::Types<D> > object
      */
      BdStepFactory(Rp::BdSimulator<D, Rpg::Types<D> >& simulator);

      /**
      * Method to create any BdStep supplied with PSCF.
      *
      * \param className name of the BdStep subclass
      * \return BdStep* pointer to new instance of className
      */
      Rp::BdStep<D, Rpg::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::BdStep<D, Rpg::Types<D> > >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Rp::BdSimulator<D, Rpg::Types<D> >* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class BdStepFactory<1>;
   extern template class BdStepFactory<2>;
   extern template class BdStepFactory<3>;

}
}
#endif
