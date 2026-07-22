/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FieldIoBase<1,CUT>;
      template class FieldIoBase<2,CUT>;
      template class FieldIoBase<3,CUT>;
      template class FieldIo<1,CUT>;
      template class FieldIo<2,CUT>;
      template class FieldIo<3,CUT>;
   }
}
