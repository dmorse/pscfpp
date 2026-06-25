/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CompressorFactory.h"  

// Subclasses of Compressor 
#include "AmCompressor.h"
#include "LrCompressor.h"
#include "LrAmCompressor.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   CompressorFactory<D>::CompressorFactory(Rp::System<D, Rpg::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Compressor subclass className.
   */
   template <int D>
   Rp::Compressor<D, Rpg::Types<D> >* 
   CompressorFactory<D>::factory(std::string const &className) const
   {
      Rp::Compressor<D, Rpg::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "Compressor" || className == "AmCompressor") {
         ptr = new AmCompressor<D>(*sysPtr_);
      } else if (className == "LrCompressor") {
         ptr = new LrCompressor<D>(*sysPtr_);
      } else if (className == "LrAmCompressor") {
         ptr = new LrAmCompressor<D>(*sysPtr_);
      }
      
      return ptr;
   }


   // Explicit instantiation definitions
   template class CompressorFactory<1>;
   template class CompressorFactory<2>;
   template class CompressorFactory<3>;

}
}
