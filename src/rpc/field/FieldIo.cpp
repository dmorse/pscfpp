/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"

namespace Pscf {
   namespace Rp {
      template class Rp::FieldIoBase<1, Rpc::Types<1> >;
      template class Rp::FieldIoBase<2, Rpc::Types<2> >;
      template class Rp::FieldIoBase<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class FieldIo<1>;
      template class FieldIo<2>;
      template class FieldIo<3>;
   }
}
