#ifndef RPC_TRAJECTORY_READER_FACTORY_TPP
#define RPC_TRAJECTORY_READER_FACTORY_TPP

#include "TrajectoryReaderFactory.h"

#include <rpc/System.h>

// Subclasses of ConfigIo
#include "RGridTrajectoryReader.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   TrajectoryReaderFactory<D>::TrajectoryReaderFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Trajectory subclass className.
   */
   template <int D>
   TrajectoryReader<D>* 
   TrajectoryReaderFactory<D>::factory(const std::string &className) const
   {
      TrajectoryReader<D> *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "RGridTrajectoryReader" 
          || className == "TrajectoryReader") {
         ptr = new RGridTrajectoryReader<D>(*sysPtr_);
      }
 
      return ptr;
   }

}
}
#endif
