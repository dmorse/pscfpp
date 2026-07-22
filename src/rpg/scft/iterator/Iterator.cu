/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>

#include <rp/scft/iterator/Iterator.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Iterator<1,CUT>;
      template class Iterator<2,CUT>;
      template class Iterator<3,CUT>;
   }
}
