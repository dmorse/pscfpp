/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CField.tpp"

namespace Pscf {
namespace Prdc {

   template class CField<1, CudaTp<1> >;
   template class CField<2, CudaTp<2> >;
   template class CField<3, CudaTp<3> >;

} // namespace Prdc
} // namespace Pscf
