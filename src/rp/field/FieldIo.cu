/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/fieldIo/fieldCheck.h>
#include <pscf/backend/cuda/VecOp.h>
#include <pscf/backend/cuda/complex.h>

#include <rp/field/FieldIo.tpp>   // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FieldIo<1,CUT>;
      template class FieldIo<2,CUT>;
      template class FieldIo<3,CUT>;
   }
}
