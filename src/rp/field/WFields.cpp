/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WFields.h"                // class template header

#include <pscf/backends/CPT.h>
#include <rpc/field/FieldIo.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/VecOp.h>

#include <rp/field/WFieldsBase.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class WFieldsBase<1,CPT>;
      template class WFieldsBase<2,CPT>;
      template class WFieldsBase<3,CPT>;
      template class WFields<1,CPT>;
      template class WFields<2,CPT>;
      template class WFields<3,CPT>;
   }
}
