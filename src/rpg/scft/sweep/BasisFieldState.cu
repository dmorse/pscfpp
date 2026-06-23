/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>

#include <rp/scft/sweep/BasisFieldState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BasisFieldState< 1, Rpg::Types<1> >;
      template class BasisFieldState< 2, Rpg::Types<2> >;
      template class BasisFieldState< 3, Rpg::Types<3> >;
   }
}
