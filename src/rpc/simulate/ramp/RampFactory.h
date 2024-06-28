#ifndef RPC_RAMP_FACTORY_H
#define RPC_RAMP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/simulate/ramp/Ramp.h>
#include <util/param/Factory.h>  
#include <string>

namespace Pscf {
namespace Rpc {

   template <int D> class Simulator;

   using namespace Util;

   /**
   * Factory for subclasses of Ramp.
   *
   * \ingroup Rpc_Simulate_Ramp_Module
   */
   template <int D>
   class RampFactory : public Factory< Ramp<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D> object
      */
      RampFactory(Simulator<D>& simulator);

      /**
      * Method to create any Ramp supplied with PSCF.
      *
      * \param className name of the Ramp subclass
      * \return Ramp* pointer to new instance of className
      */
      Ramp<D>* factory(const std::string & className) const;

      using Factory< Ramp<D> >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Simulator<D>* simulatorPtr_;

   };

   #ifndef RPC_RAMP_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class RampFactory<1>;
   extern template class RampFactory<2>;
   extern template class RampFactory<3>;
   #endif

}
}
#endif
