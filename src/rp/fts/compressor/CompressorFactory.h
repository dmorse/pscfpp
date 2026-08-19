#ifndef RP_COMPRESSOR_FACTORY_H
#define RP_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>             // base class template
#include <rp/fts/compressor/Compressor.h>   // base class argument
#include <pscf/backends/TmplDeclare.h>      // preprocessor macros

#include <string>

namespace Pscf {
namespace Rp {

   // Forward declaration
   template <int D, class T> class System;

   using namespace Util;

   /**
   * Factory for subclasses of Compressor.
   *
   * \ingroup Rp_Fts_Compressor_Module
   */

   template <int D, class T>
   class CompressorFactory 
     : public Factory< Compressor<D,T> > 
   {

   public:

      /// Constructor
      CompressorFactory(System<D,T>& system);

      /**
      * Method to create any Compressor supplied with PSCF.
      *
      * \param className name of the Compressor subclass
      * \return Compressor* pointer to new instance of className
      */
      Compressor<D,T>* factory(const std::string &className) const;

      using Factory< Compressor<D,T> >::trySubfactories;

   private:

      /// Pointer to the parent system.
      System<D,T>* sysPtr_;

      /// Alias to base class
      using FactoryT = Factory< Compressor<D,T> >;

   };

   // Explicit instantation declarations
   PSCF_TMPL_DECLARE(CompressorFactory)

}
}
#endif
