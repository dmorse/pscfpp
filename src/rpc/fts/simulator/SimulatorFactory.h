#ifndef RPC_SIMULATOR_FACTORY_H
#define RPC_SIMULATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpc/fts/simulator/Simulator.h>

#include <string>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Factory for subclasses of Simulator.
   *
   * \ingroup Rpc_Fts_Simulator_Module
   */
   template <int D>
   class SimulatorFactory : public Factory< Rp::Simulator<D, Rpc::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent Rp::System<D, Rpc::Types<D> > object
      */
      SimulatorFactory(Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Method to create any Simulator supplied with PSCF.
      *
      * \param className name of the Simulator subclass
      * \return Simulator* pointer to new instance of className
      */
      Rp::Simulator<D, Rpc::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::Simulator<D, Rpc::Types<D> > >::trySubfactories;
      using Factory< Rp::Simulator<D, Rpc::Types<D> > >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      Rp::System<D, Rpc::Types<D> >* systemPtr_;

   };

   // Explicit instantiation declarations
   extern template class SimulatorFactory<1>;
   extern template class SimulatorFactory<2>;
   extern template class SimulatorFactory<3>;

}
}
#endif
