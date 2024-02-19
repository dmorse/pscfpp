#ifndef R1D_ITERATOR_FACTORY_H
#define R1D_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <r1d/iterator/Iterator.h>
#include <r1d/System.h>

#include <string>

namespace Pscf {
namespace R1d {

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * The default iterator, invoked by className Iterator, is the
   * Anderson mixing iterator (className AmIterator).
   * 
   * \ingroup R1d_Iterator_Module
   */
   class IteratorFactory : public Factory<Iterator> 
   {

   public:

      /**
      * Constructor
      *
      * \param system  parent System object
      */
      IteratorFactory(System& system);

      /**
      * Method to create any Iterator supplied with PSCF.
      *
      * \param className name of the Iterator subclass
      * \return Iterator* pointer to new instance of className
      */
      Iterator* factory(const std::string &className) const;

      using Factory< Iterator >::trySubfactories;

   private:

      /// Pointer to the parent system.
      System* sysPtr_;

   };


}
}
#endif
