#ifndef RPC_BD_STEP_FACTORY_H
#define RPC_BD_STEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpc/fts/brownian/BdStep.h>

#include <string>

namespace Pscf {
namespace Rpc {

   // Forward declaration
   template <int D> class BdSimulator;

   using namespace Util;

   /**
   * Factory for subclasses of BdStep.
   *
   * \ingroup Rpc_Fts_Brownian_Module
   */
   template <int D>
   class BdStepFactory : public Factory< Rp::BdStep<D, Rpc::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rp::BdSimulator<D, Rpc::Types<D> > object
      */
      BdStepFactory(Rp::BdSimulator<D, Rpc::Types<D> >& simulator);

      /**
      * Method to create any BdStep supplied with PSCF.
      *
      * \param className name of the BdStep subclass
      * \return BdStep* pointer to new instance of className
      */
      Rp::BdStep<D, Rpc::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::BdStep<D, Rpc::Types<D> > >::trySubfactories;

   private:

      /// Pointer to the parent simulator.
      Rp::BdSimulator<D, Rpc::Types<D> >* simulatorPtr_;

   };

   // Explicit instantiation declarations
   extern template class BdStepFactory<1>;
   extern template class BdStepFactory<2>;
   extern template class BdStepFactory<3>;

}
}
#endif
