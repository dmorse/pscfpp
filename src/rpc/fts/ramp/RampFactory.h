#ifndef RPC_RAMP_FACTORY_H
#define RPC_RAMP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/RampFactory.h>
#include <rpc/system/Types.h>

#if 0
#include <rpc/fts/ramp/Ramp.h>
#include <util/param/Factory.h>  
#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Factory for subclasses of Ramp.
   *
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class RampFactory : public Factory< Rp::Ramp<D, Rpc::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rp::Simulator<D, Rpc::Types<D> > object
      */
      RampFactory(Rp::Simulator<D, Rpc::Types<D> >& simulator);

      /**
      * Method to create any Ramp supplied with PSCF.
      *
      * \param className name of the Ramp subclass
      * \return Ramp* pointer to new instance of className
      */
      Rp::Ramp<D, Rpc::Types<D> >* factory(const std::string & className) const;

      using Factory< Rp::Ramp<D, Rpc::Types<D> > >::trySubfactories;

   private:
      
      /// Pointer to the parent simulator.
      Rp::Simulator<D, Rpc::Types<D> >* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class Rp::RampFactory<1, Rpc::Types<1> >;
   extern template class Rp::RampFactory<2, Rpc::Types<2> >;
   extern template class Rp::RampFactory<3, Rpc::Types<3> >;

}
}
#endif

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RampFactory<1, Rpc::Types<1> >;
      extern template class RampFactory<2, Rpc::Types<2> >;
      extern template class RampFactory<3, Rpc::Types<3> >;
   }
}
#endif
