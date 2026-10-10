/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/fieldIo/fieldCheck.h>
#include <pscf/backend/cpp/HostRecvArray.h>
#include <pscf/backend/cpp/HostArray.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/complex.h>

#include <rp/field/FieldIo.tpp>   // base class implementation

// Explicit specialization definitions
namespace Pscf {
   namespace Rp {
      template class Rp::FieldIo<1,CPT>;
      template class Rp::FieldIo<2,CPT>;
      template class Rp::FieldIo<3,CPT>;
   }
}
