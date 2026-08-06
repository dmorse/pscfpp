/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>

#include <rp/scft/sweep/BasisFieldState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BasisFieldState<1,CPT>;
      template class BasisFieldState<2,CPT>;
      template class BasisFieldState<3,CPT>;
   }
}
