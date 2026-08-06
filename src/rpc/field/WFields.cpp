/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/FieldIo.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/VecOp.h>

#include <rp/field/WFieldsBase.tpp> // base class implementation
#include <rp/field/WFields.h>       // class header

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
