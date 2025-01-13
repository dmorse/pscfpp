/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"
#include <util/math/Constants.h>

namespace Pscf {

namespace Prdc {

   using namespace Cpu;

   // Explicit instantiation
   template class FieldIoReal<1, RField<1>, RFieldDft<1>, FFT<1> >;
   template class FieldIoReal<2, RField<2>, RFieldDft<2>, FFT<2> >;
   template class FieldIoReal<3, RField<3>, RFieldDft<3>, FFT<3> >;
}

namespace Rpc {

   // Explicit instantiations
   template class FieldIo<1>;
   template class FieldIo<2>;
   template class FieldIo<3>;

} 

} // namespace Pscf
