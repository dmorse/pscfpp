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
   class SweepFactory : public Factory< Sweep<D> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      SweepFactory(System<D>& system);

      /**
      * Method to create any Sweep subclass.
      *
      * \param className name of the Sweep subclass
      * \return Sweep<D>* pointer to new instance of speciesName
      */
      Sweep<D>* factory(std::string const & className) const;

      using Factory< Sweep<D> >::trySubfactories;
      using Factory< Sweep<D> >::readObjectOptional;

   private:

      // Pointer to parent system object.
      System<D>* systemPtr_;

   };

   // Explicit instantiation declarations
   extern template class SweepFactory<1>;
   extern template class SweepFactory<2>;
   extern template class SweepFactory<3>;

}
}
#endif
