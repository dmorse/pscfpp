#ifndef RPG_COMPRESSOR_FACTORY_H
#define RPG_COMPRESSOR_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>            // base class template
#include <rp/fts/compressor/Compressor.h>  // base template argument
#include <rpg/system/Types.h>              // argument of argument

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
   }
}

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
   class CompressorFactory : public Factory< Rp::Compressor<D, Rpg::Types<D> > > 
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent system
      */
      CompressorFactory(Rp::System<D, Rpg::Types<D> >& system);

      /**
      * Method to create a Compressor.
      *
      * \param className name of the Compressor subclass
      * \return Compressor* pointer to new instance of className
      */
      Rp::Compressor<D, Rpg::Types<D> >* factory(std::string const & className) 
      const;

      using Factory< Rp::Compressor<D, Rpg::Types<D> > >::trySubfactories;

   private:

      /// Pointer to the parent system.
      Rp::System<D, Rpg::Types<D> >* sysPtr_;

   };

   // Explicit instantiation declarations
   extern template class CompressorFactory<1>;
   extern template class CompressorFactory<2>;
   extern template class CompressorFactory<3>;

}
}
#endif
