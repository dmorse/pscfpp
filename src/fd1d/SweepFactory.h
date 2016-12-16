#ifndef FD1D_SWEEP_FACTORY_H
#define FD1D_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include "Sweep.h"

#include <string>

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   /**
   * Default Factory for subclasses of Sweep.
   */
   class SweepFactory : public Factory<Sweep> 
   {

   public:

      /**
      * Method to create any species supplied with Simpatico.
      *
      * \param speciesName name of the Sweep subclass
      * \return Sweep* pointer to new instance of speciesName
      */
      Sweep* factory(const std::string &speciesName) const;

   };

}
}
#endif
