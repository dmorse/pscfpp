/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/WFields.h>

#include <rp/scft/sweep/BasisFieldState.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BasisFieldState< 1, CppTp<1> >;
      template class BasisFieldState< 2, CppTp<2> >;
      template class BasisFieldState< 3, CppTp<3> >;
   }
}
