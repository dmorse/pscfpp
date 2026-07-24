/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/CFields.h>
#include <rp/field/FieldIo.h>
#include <prdc/field/cpu/RField.h>
#include <rp/field/CFields.tpp>   // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class CFields<1,CPT>;
      template class CFields<2,CPT>;
      template class CFields<3,CPT>;
   }
}
