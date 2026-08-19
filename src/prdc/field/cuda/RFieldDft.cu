/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.tpp"

// Explicit instantiation definitions
namespace Pscf {
namespace Prdc {
   template class RFieldDft<1,CUT>;
   template class RFieldDft<2,CUT>;
   template class RFieldDft<3,CUT>;
}
}
