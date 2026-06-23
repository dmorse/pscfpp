#ifndef RPG_ITERATOR_FACTORY_H
#define RPG_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/iterator/Iterator.h>
#include <rpg/system/Types.h>
#include <util/param/Factory.h>  

#include <string>

namespace Pscf {
namespace Rp {

   // Forward declaration
   template <int D, class T> class IteratorFactory; 

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
      IteratorFactory(Rp::System<D, Rpg::Types<D> >& system);

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
      Rp::System<D, Rpg::Types<D> >* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class IteratorFactory<1, Rpg::Types<1> >;
   extern template class IteratorFactory<2, Rpg::Types<2> >;
   extern template class IteratorFactory<3, Rpg::Types<3> >;

}
}
#endif
