#ifndef PSPG_ITERATOR_FACTORY_H
#define PSPG_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/iterator/Iterator.h>
#include <rpg/System.h>

#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rpg_Iterator_Module
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

   private:

      /// Pointer to the system object.
      System<D>* sysPtr_;

   };

   #ifndef PSPG_ITERATOR_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class IteratorFactory<1>;
   extern template class IteratorFactory<2>;
   extern template class IteratorFactory<3>;
   #endif

}
}
#endif
