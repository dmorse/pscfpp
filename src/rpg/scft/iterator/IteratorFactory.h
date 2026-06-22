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
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */

   template <int D>
   class IteratorFactory : public Factory< Rp::Iterator<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system
      */
      IteratorFactory(System<D>& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Rp::Iterator<D, Rpg::Types<D> >* factory(const std::string &className) const;

      // Inherited member functions
      using Factory< Rp::Iterator<D, Rpg::Types<D> > >::trySubfactories;
      using Factory< Rp::Iterator<D, Rpg::Types<D> > >::readObjectOptional;

   private:

      /// Pointer to the parent system object.
      System<D>* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class IteratorFactory<1>;
   extern template class IteratorFactory<2>;
   extern template class IteratorFactory<3>;

}
}
#endif
