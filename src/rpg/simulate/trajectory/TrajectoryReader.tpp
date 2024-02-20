#ifndef PSPG_TRAJECTORY_READER_TPP
#define PSPG_TRAJECTORY_READER_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"

namespace Pscf {
namespace Rpg 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   TrajectoryReader<D>::TrajectoryReader(System<D>& system)
    : systemPtr_(&system)
   {}
   

  
}
}
#endif
