/*
* PSCF - Polymer Self-Consistent HostField Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/HostField.tpp>

namespace Pscf {
namespace Prdc { 
namespace Cuda {

   template class HostField<cudaReal>;
   template class HostField<cudaComplex>;

}
}
}
