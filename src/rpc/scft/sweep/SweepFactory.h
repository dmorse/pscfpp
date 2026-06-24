#ifndef RPC_SWEEP_FACTORY_H
#define RPC_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include "Sweep.h"

#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   
   /**
   * Default Factory for subclasses of Sweep.
   *
   * \ingroup Rpc_Scft_Sweep_Module
   */
   template <int D>
   class SweepFactory : public Factory< Rp::Sweep<D, Rpc::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      SweepFactory(Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Method to create any Sweep subclass.
      *
      * \param className name of the Sweep subclass
      * \return Rp::Sweep<D, Rpc::Types<D> >* pointer to new instance of speciesName
      */
      Rp::Sweep<D, Rpc::Types<D> >* factory(std::string const & className) const;

      using Factory< Rp::Sweep<D, Rpc::Types<D> > >::trySubfactories;
      using Factory< Rp::Sweep<D, Rpc::Types<D> > >::readObjectOptional;

   private:

      // Pointer to parent system object.
      Rp::System<D, Rpc::Types<D> >* systemPtr_;

   };

   // Explicit instantiation declarations
   extern template class SweepFactory<1>;
   extern template class SweepFactory<2>;
   extern template class SweepFactory<3>;

}
}
#endif
