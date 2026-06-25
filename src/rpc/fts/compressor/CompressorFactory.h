#ifndef RPC_COMPRESSOR_FACTORY_H
#define RPC_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>            // base class template
#include <rp/fts/compressor/Compressor.h>  // base template argument
#include <rpc/system/Types.h>              // argument of argument

#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Factory for subclasses of Compressor.
   *
   * \ingroup Rpc_Fts_Compressor_Module
   */

   template <int D>
   class CompressorFactory : public Factory< Rp::Compressor<D, Rpc::Types<D> > > 
   {

   public:

      /// Constructor
      CompressorFactory(Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Method to create any Compressor supplied with PSCF.
      *
      * \param className name of the Compressor subclass
      * \return Compressor* pointer to new instance of className
      */
      Rp::Compressor<D, Rpc::Types<D> >* factory(const std::string &className) const;

      using Factory< Rp::Compressor<D, Rpc::Types<D> > >::trySubfactories;

   private:

      /// Pointer to the parent system.
      Rp::System<D, Rpc::Types<D> >* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class CompressorFactory<1>;
   extern template class CompressorFactory<2>;
   extern template class CompressorFactory<3>;

}
}
#endif
