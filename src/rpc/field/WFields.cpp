/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                // class header
#include <rpc/field/FieldIo.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <rp/field/WFieldsBase.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class WFieldsBase<1, Cpp<1> >;
      template class WFieldsBase<2, Cpp<2> >;
      template class WFieldsBase<3, Cpp<3> >;
      template class WFields<1, Cpp<1> >;
      template class WFields<2, Cpp<2> >;
      template class WFields<3, Cpp<3> >;
   }
}
