/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"
#include <rpc/field/FieldIo.h>
#include <prdc/field/cpu/RField.h>
#include <rp/field/CFields.tpp>   // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      template class CFields<1, Cpp<1> >;
      template class CFields<2, Cpp<2> >;
      template class CFields<3, Cpp<3> >;
   }
}
