/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mask.tpp"
#include <prdc/field/MaskTmpl.tpp> // needed for implementation

namespace Pscf {

namespace Prdc {

   // Explicit instantiation
   template class MaskTmpl< 1, Cuda::RField<1>, Rpg::FieldIo<1> >;
   template class MaskTmpl< 2, Cuda::RField<2>, Rpg::FieldIo<2> >;
   template class MaskTmpl< 3, Cuda::RField<3>, Rpg::FieldIo<3> >;

} // namespace Prdc

namespace Rpg {

   template class Mask<1>;
   template class Mask<2>;
   template class Mask<3>;
   
} // namespace Rpg

} // namespace Pscf
