/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.tpp"
#include <prdc/field/MaskTmpl.tpp> // needed for implementation

namespace Pscf {

namespace Prdc {

   // Explicit instantiation
   template class MaskTmpl< 1, Rpc::FieldIo<1>, Cpu::RField<1> >;
   template class MaskTmpl< 2, Rpc::FieldIo<2>, Cpu::RField<2> >;
   template class MaskTmpl< 3, Rpc::FieldIo<3>, Cpu::RField<3> >;

} // namespace Prdc

namespace Rpc {

   template class Mask<1>;
   template class Mask<2>;
   template class Mask<3>;
   
} // namespace Rpc

} // namespace Pscf
