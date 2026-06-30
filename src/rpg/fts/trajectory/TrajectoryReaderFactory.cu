/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.h"

// Subclasses of ConfigIo
#include "RGridTrajectoryReader.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   TrajectoryReaderFactory<D>::TrajectoryReaderFactory(Rp::System<D, Rpg::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of TrajectoryReader subclass className.
   */
   template <int D>
   Rp::TrajectoryReader<D, Rpg::Types<D> >* 
   TrajectoryReaderFactory<D>::factory(const std::string &className) const
   {
      Rp::TrajectoryReader<D, Rpg::Types<D> > *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "RGridTrajectoryReader" 
          || className == "TrajectoryReader") {
        ptr = new Rp::RGridTrajectoryReader<D, Rpg::Types<D> >(*sysPtr_);
      } 
      return ptr;
   }

   // Explicit instantiation declarations
   template class TrajectoryReaderFactory<1>;
   template class TrajectoryReaderFactory<2>;
   template class TrajectoryReaderFactory<3>;

}
}

