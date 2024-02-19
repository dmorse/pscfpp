#ifndef R1D_SWEEP_FACTORY_H
#define R1D_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include "Sweep.h"

#include <string>

namespace Pscf {
namespace R1d {

   using namespace Util;

   /**
   * Default Factory for subclasses of Sweep.
   *
   * \ingroup R1d_Sweep_Module
   */
   class SweepFactory : public Factory<Sweep> 
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      SweepFactory(System& system);

      /**
      * Method to create any Sweep subclass.
      *
      * \param className name of the Sweep subclass
      * \return Sweep* pointer to new instance of speciesName
      */
      Sweep* factory(std::string const & className) const;

   private:

      System* systemPtr_;

   };

}
}
#endif
