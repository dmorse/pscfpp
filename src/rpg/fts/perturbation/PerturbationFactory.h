#ifndef RPG_PERTURBATION_FACTORY_H
#define RPG_PERTURBATION_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/fts/perturbation/Perturbation.h>
#include <util/param/Factory.h>  
#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Factory for subclasses of Perturbation.
   *
   * \ingroup Rpg_Fts_Perturbation_Module
   */
   template <int D>
   class PerturbationFactory : public Factory< Rp::Perturbation<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rp::Simulator<D, Rpg::Types<D> > object
      */
      PerturbationFactory(Rp::Simulator<D, Rpg::Types<D> >& simulator);

      /**
      * Method to create any Perturbation supplied with PSCF.
      *
      * \param className name of the Perturbation subclass
      * \return Perturbation* pointer to new instance of className
      */
      Rp::Perturbation<D, Rpg::Types<D> >* factory(const std::string & className) const;

      using Factory< Rp::Perturbation<D, Rpg::Types<D> > >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Rp::Simulator<D, Rpg::Types<D> >* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class PerturbationFactory<1>;
   extern template class PerturbationFactory<2>;
   extern template class PerturbationFactory<3>;

}
}
#endif
