#ifndef RPC_PERTURBATION_FACTORY_H
#define RPC_PERTURBATION_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/perturbation/Perturbation.h>
#include <util/param/Factory.h>  
#include <string>

namespace Pscf {
namespace Rpc {

   template <int D> class Simulator;

   using namespace Util;

   /**
   * Factory for subclasses of Perturbation.
   *
   * \ingroup Rpc_Fts_Perturbation_Module
   */
   template <int D>
   class PerturbationFactory : public Factory< Perturbation<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D> object
      */
      PerturbationFactory(Simulator<D>& simulator);

      /**
      * Method to create any Perturbation supplied with PSCF.
      *
      * \param className name of the Perturbation subclass
      * \return Perturbation* pointer to new instance of className
      */
      Perturbation<D>* factory(const std::string & className) const;

      using Factory< Perturbation<D> >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Simulator<D>* simulatorPtr_;

   };

   #ifndef RPC_PERTURBATION_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class PerturbationFactory<1>;
   extern template class PerturbationFactory<2>;
   extern template class PerturbationFactory<3>;
   #endif

}
}
#endif
