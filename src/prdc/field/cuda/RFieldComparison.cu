/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.tpp"

namespace Pscf {
namespace Prdc {

   template class RFieldComparison<1, CudaTp<1> >;
   template class RFieldComparison<2, CudaTp<2> >;
   template class RFieldComparison<3, CudaTp<3> >;

} // namespace Pscf::Prdc
} // namespace Pscf
