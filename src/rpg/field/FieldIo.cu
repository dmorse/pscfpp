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
      template class FieldIoBase<1, CudaTp<1> >;
      template class FieldIoBase<2, CudaTp<2> >;
      template class FieldIoBase<3, CudaTp<3> >;
      template class FieldIo<1, CudaTp<1> >;
      template class FieldIo<2, CudaTp<2> >;
      template class FieldIo<3, CudaTp<3> >;
   }
}
