/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.tpp"

// Explicit instantiation definitions
namespace Pscf {
namespace Prdc {
   template class RField<1, CudaTp<1> >;
   template class RField<2, CudaTp<2> >;
   template class RField<3, CudaTp<3> >;
}
}
