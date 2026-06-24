/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BdStep.h"
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/system/System.h>

#include <rp/fts/brownian/BdStep.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdStep< 1, Rpg::Types<1> >;
      template class BdStep< 2, Rpg::Types<2> >;
      template class BdStep< 3, Rpg::Types<3> >;
   }
}
