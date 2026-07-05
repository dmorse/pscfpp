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
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   CompressorFactory<D>::CompressorFactory(Rp::System<D, Rpc::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Compressor subclass className.
   */
   template <int D>
   Rp::Compressor<D, Rpc::Types<D> >* 
   CompressorFactory<D>::factory(std::string const &className) const
   {
      Rp::Compressor<D, Rpc::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "Compressor" || className == "LrAmCompressor") {
         ptr = new Rp::LrAmCompressor<D, Rpc::Types<D>, DArray<double> >(*sysPtr_);
      } else if (className == "AmCompressor") {
         ptr = new Rp::AmCompressor<D, Rpc::Types<D>, DArray<double> >(*sysPtr_);
      } else if (className == "LrCompressor") {
         ptr = new Rp::LrCompressor<D, Rpc::Types<D> >(*sysPtr_);
      }
  
      return ptr;
   }

   // Explicit instantiation definitions
   template class CompressorFactory<1>;
   template class CompressorFactory<2>;
   template class CompressorFactory<3>;

}
}
