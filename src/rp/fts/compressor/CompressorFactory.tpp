/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CompressorFactory.h"  

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   CompressorFactory<D,T>::CompressorFactory(System<D,T>& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Compressor subclass className.
   */
   template <int D, class T>
   Compressor<D,T>* 
   CompressorFactory<D,T>::factory(std::string const &className) const
   {
      Compressor<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "Compressor" || className == "LrAmCompressor") {
         ptr = new LrAmCompressor<D,T>(*sysPtr_);
      } else if (className == "AmCompressor") {
         ptr = new AmCompressor<D,T>(*sysPtr_);
      } else if (className == "LrCompressor") {
         ptr = new LrCompressor<D,T>(*sysPtr_);
      }
  
      return ptr;
   }

}
}
