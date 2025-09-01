#ifndef RPG_COMPRESSOR_FACTORY_H
#define RPG_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/system/System.h>

#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Factory for subclasses of Compressor.
   *
   * \ingroup Rpg_Fts_Compressor_Module
   */

   template <int D>
   class CompressorFactory : public Factory< Compressor<D> > 
   {

   public:

      /// Constructor
      CompressorFactory(System<D>& system);

      /**
      * Method to create any Compressor supplied with PSCF.
      *
      * \param className name of the Compressor subclass
      * \return Compressor* pointer to new instance of className
      */
      Compressor<D>* factory(const std::string &className) const;

      using Factory< Compressor<D> >::trySubfactories;

   private:

      /// Pointer to the parent system.
      System<D>* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class CompressorFactory<1>;
   extern template class CompressorFactory<2>;
   extern template class CompressorFactory<3>;

}
}
#endif
