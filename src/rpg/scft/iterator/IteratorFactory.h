#ifndef RPG_ITERATOR_FACTORY_H
#define RPG_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/scft/iterator/Iterator.h>
#include <rpg/system/System.h>
#include <util/param/Factory.h>  

#include <string>

namespace Pscf {
   namespace Rp {
      template <int D, class T> class IteratorFactory;
   }
   namespace Rpg {
      template <int D> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */

   template <int D>
   class IteratorFactory<D, Rpg::Types<D> > 
     : public Factory< Iterator<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system
      */
      IteratorFactory(Rpg::System<D>& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator<D, Rpg::Types<D> >* factory(const std::string &className) const;

      // Inherited member functions
      using Factory< Iterator<D, Rpg::Types<D> > >::trySubfactories;
      using Factory< Iterator<D, Rpg::Types<D> > >::readObjectOptional;

   private:

      /// Pointer to the parent system object.
      Rpg::System<D>* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class IteratorFactory<1, Rpg::Types<1> >;
   extern template class IteratorFactory<2, Rpg::Types<2> >;
   extern template class IteratorFactory<3, Rpg::Types<3> >;

}
}
#endif
