/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"

// Subclasses of Sweep
#include "LinearSweep.h"

namespace Pscf{
namespace Rpc{

   using namespace Util;

   template <int D>
   SweepFactory<D>::SweepFactory(Rp::System<D, Rpc::Types<D> >& system)
    : systemPtr_(&system)
   {}

   /**
   * Return a pointer to a Sweep subclass with name className
   */
   template <int D>
   Sweep<D>* SweepFactory<D>::factory(std::string const & className) 
   const
   {
       Sweep<D> *ptr = nullptr;

       // First check if name is known by any subfactories
       ptr = trySubfactories(className);
       if (ptr) return ptr;

       // Explicit class names
       if (className == "Sweep" || className == "LinearSweep") {
           ptr = new LinearSweep<D>(*systemPtr_);
       }

       return ptr;
   }

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rpc {
      template class SweepFactory<1>;
      template class SweepFactory<2>;
      template class SweepFactory<3>;
   } 
} 
