/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IteratorMediatorCUDA.tpp"

namespace Pscf {
namespace Pspg {
   template class IteratorMediatorCUDA<1>;
   template class IteratorMediatorCUDA<2>;
   template class IteratorMediatorCUDA<3>;
}
}
