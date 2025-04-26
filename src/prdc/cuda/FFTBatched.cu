/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.tpp"

namespace Pscf {
namespace Prdc {
namespace Cuda {

   template class FFTBatched<1>;
   template class FFTBatched<2>;
   template class FFTBatched<3>;

} // namespace Pscf::Prdc::Cuda
} // namespace Pscf::Prdc
} // namespace Pscf
