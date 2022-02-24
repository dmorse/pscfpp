#ifndef PSPC_ITERATOR_FACTORY_H
#define PSPC_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <pscf/iterator/Iterator.h>
#include <pspc/iterator/IteratorMediatorCPU.h>

#include <string>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */

   template <int D>
   class IteratorFactory : public Factory< Iterator<FieldCPU> > 
   {

   public:

      /// Constructor
      IteratorFactory(IteratorMediatorCPU<D>& iterMed);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator<FieldCPU>* factory(const std::string &className) const;

      using Factory< Iterator<FieldCPU> >::trySubfactories;

   private:

      /// Pointer to an IteratorMediatorCPU object.
      IteratorMediatorCPU<D>* iterMedPtr_;

   };

   #ifndef PSPC_ITERATOR_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class IteratorFactory<1>;
   extern template class IteratorFactory<2>;
   extern template class IteratorFactory<3>;
   #endif

}
}
#endif
