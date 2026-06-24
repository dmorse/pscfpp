/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"

#include "LinearSweep.h"

namespace Pscf{
namespace Rpg{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   SweepFactory<D>::SweepFactory(Rp::System<D, Rpg::Types<D> >& system)
    : systemPtr_(&system)
   {}

   /*
   * Return a pointer to a Sweep subclass with name className
   */
   template <int D>
   Rp::Sweep<D, Rpg::Types<D> >* SweepFactory<D>::factory(std::string const & className) const
   {
       Rp::Sweep<D, Rpg::Types<D> > *ptr = 0;

       // First check if name is known by any subfactories
       ptr = trySubfactories(className);
       if (ptr) return ptr;

       // Explicit class names
       if (className == "Sweep" || className == "LinearSweep") {
           ptr = new Rp::LinearSweep<D, Rpg::Types<D> >(*systemPtr_);
       }

       return ptr;
   }

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Rpg {
      template class SweepFactory<1>;
      template class SweepFactory<2>;
      template class SweepFactory<3>;
   }
}
