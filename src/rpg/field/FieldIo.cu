/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"
#include <util/math/Constants.h>

namespace Pscf {

namespace Prdc {

   using namespace Pscf::Prdc::Cuda;

   // Explicit instantiation
   template class FieldIoReal<1, Prdc::Cuda::RField<1>, Prdc::Cuda::RFieldDft<1>, Prdc::Cuda::FFT<1> >;
   template class FieldIoReal<2, Prdc::Cuda::RField<2>, Prdc::Cuda::RFieldDft<2>, Prdc::Cuda::FFT<2> >;
   template class FieldIoReal<3, Prdc::Cuda::RField<3>, Prdc::Cuda::RFieldDft<3>, Prdc::Cuda::FFT<3> >;

}

namespace Rpg {

   // Explicit instantiations
   template class FieldIo<1>;
   template class FieldIo<2>;
   template class FieldIo<3>;

} 

} // namespace Pscf
