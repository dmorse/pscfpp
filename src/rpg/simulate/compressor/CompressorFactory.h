#ifndef RPG_COMPRESSOR_FACTORY_H
#define RPG_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <rpg/simulate/compressor/Compressor.h>
#include <rpg/System.h>

#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Factory for subclasses of Compressor.
   *
   * \ingroup Rpg_Compressor_Module
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

   #ifndef RPG_COMPRESSOR_FACTORY_TPP
   // Suppress implicit instantiation
   extern template class CompressorFactory<1>;
   extern template class CompressorFactory<2>;
   extern template class CompressorFactory<3>;
   #endif

}
}
#endif
