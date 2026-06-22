#ifndef RPC_ITERATOR_FACTORY_H
#define RPC_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/scft/iterator/Iterator.h>
#include <rpc/system/Types.h>
#include <util/param/Factory.h>  

#include <string>

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class IteratorFactory;
   }
   namespace Rpc {
      template <int D> class System;
   }
}

namespace Pscf {
namespace Rp {


   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rpc_Scft_Iterator_Module
   */

   template <int D>
   class IteratorFactory<D, Rpc::Types<D> > 
    : public Factory< Rp::Iterator<D, Rpc::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      IteratorFactory(Rpc::System<D>& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className  name of the Rp::Iterator<D, Rpc::Types<D> > subclass
      * \return  pointer to new instance of className
      */
      Rp::Iterator<D, Rpc::Types<D> >* factory(const std::string &className) const;

      // Inherited member functions
      using Factory< Rp::Iterator<D, Rpc::Types<D> > >::trySubfactories;
      using Factory< Rp::Iterator<D, Rpc::Types<D> > >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      Rpc::System<D>* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class IteratorFactory<1, Rpc::Types<1> >;
   extern template class IteratorFactory<2, Rpc::Types<2> >;
   extern template class IteratorFactory<3, Rpc::Types<3> >;

}
}
#endif
