/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.h"
//#include <rpc/system/System.h>

// Subclasses of TrajectoryReader
#include "RGridTrajectoryReader.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   TrajectoryReaderFactory<D>::TrajectoryReaderFactory(Rp::System<D, Rpc::Types<D> >& system)
    : sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Trajectory subclass className.
   */
   template <int D>
   Rp::TrajectoryReader<D, Rpc::Types<D> >* 
   TrajectoryReaderFactory<D>::factory(const std::string &className) const
   {
      Rp::TrajectoryReader<D, Rpc::Types<D> > *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "RGridTrajectoryReader" 
          || className == "TrajectoryReader") {
         ptr = new Rp::RGridTrajectoryReader<D, Rpc::Types<D> >(*sysPtr_);
      }
 
      return ptr;
   }

   // Explicit instantiation definitions
   template class TrajectoryReaderFactory<1>;
   template class TrajectoryReaderFactory<2>;
   template class TrajectoryReaderFactory<3>;

}
}
