#ifndef PSPG_ITERATOR_FACTORY_H
#define PSPG_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <pscf/iterator/Iterator.h>
#include <pspg/iterator/IteratorMediatorCUDA.h>

#include <string>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Pspg_Iterator_Module
   */

   template <int D>
   class IteratorFactory : public Factory< Iterator<FieldCUDA> > 
   {

   public:

      /// Constructor
      IteratorFactory(IteratorMediatorCUDA<D>& iterMed);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator<FieldCUDA>* factory(const std::string &className) const;

      using Factory< Iterator<FieldCUDA> >::trySubfactories;

   private:

      /// Pointer to an IteratorMediatorCUDA object.
      IteratorMediatorCUDA<D>* iterMedPtr_;

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
