/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.tpp"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   template class RFieldDft<1>;
   template class RFieldDft<2>;
   template class RFieldDft<3>;

} // namespace Pscf::Prdc::Cuda
} // namespace Pscf::Prdc
} // namespace Pscf
