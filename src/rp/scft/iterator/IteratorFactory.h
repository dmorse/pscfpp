#ifndef RP_ITERATOR_FACTORY_H
#define RP_ITERATOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>          // base class template
#include <rp/scft/iterator/Iterator.h>   // class template argument
#include <pscf/backends/TmplDeclare.h>   // declaration macros

#include <string>

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Factory for subclasses of Iterator.
   *
   * \ingroup Rp_Scft_Iterator_Module
   */

   template <int D, class T>
   class IteratorFactory : public Factory< Iterator<D,T> > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System
      */
      IteratorFactory(System<D,T>& system);

      /**
      * Create an instance of a specified subclass of Iterator<D,T>.
      *
      * \param className  name of the Iterator<D,T> subclass
      * \return  pointer to new instance of subclass className
      */
      Iterator<D,T>* factory(const std::string &className) const;

      // Inherited member functions
      using Factory< Iterator<D,T> >::trySubfactories;
      using Factory< Iterator<D,T> >::readObjectOptional;

   private:

      /// Pointer to the parent system.
      System<D,T>* sysPtr_;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(IteratorFactory)

}
}
#endif
