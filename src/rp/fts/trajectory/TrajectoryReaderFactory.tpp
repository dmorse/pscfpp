/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.h"

// Subclasses of TrajectoryReader
#include <rp/fts/trajectory/RGridTrajectoryReader.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   TrajectoryReaderFactory<D,T>::TrajectoryReaderFactory(
                                                 System<D,T>& system)
    : sysPtr_(&system)
   {}

   /*
   * Return pointer to an instance of TrajectoryReader subclass className.
   */
   template <int D, class T>
   TrajectoryReader<D,T>* 
   TrajectoryReaderFactory<D,T>::factory(const std::string &className) const
   {
      TrajectoryReader<D,T>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "RGridTrajectoryReader" 
          || className == "TrajectoryReader") {
         ptr = new RGridTrajectoryReader<D,T>(*sysPtr_);
      }
 
      return ptr;
   }

}
}
