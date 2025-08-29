#ifndef RPC_ITERATOR_FACTORY_H
#define RPC_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/scft/iterator/Iterator.h>
#include <util/param/Factory.h>  

#include <string>

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rpc_Scft_Iterator_Module
   */

   template <int D>
   class IteratorFactory : public Factory< Iterator<D> > 
   {

   public:

      /// Constructor
      IteratorFactory(System<D>& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator<D>* factory(const std::string &className) const;

      using Factory< Iterator<D> >::trySubfactories;
      using Factory< Iterator<D> >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      System<D>* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class IteratorFactory<1>;
   extern template class IteratorFactory<2>;
   extern template class IteratorFactory<3>;

}
}
#endif
