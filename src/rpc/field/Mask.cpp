/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.tpp"
#include <prdc/field/MaskReal.tpp> // needed for implementation

namespace Pscf {

namespace Prdc {

   // Explicit instantiation
   template class MaskReal< 1, Rpc::FieldIo<1>, Cpu::RField<1> >;
   template class MaskReal< 2, Rpc::FieldIo<2>, Cpu::RField<2> >;
   template class MaskReal< 3, Rpc::FieldIo<3>, Cpu::RField<3> >;

} // namespace Prdc

namespace Rpc {

   template class Mask<1>;
   template class Mask<2>;
   template class Mask<3>;
   
} // namespace Rpc

} // namespace Pscf
