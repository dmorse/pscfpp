/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDftComparison.tpp"

namespace Pscf {
   namespace Prdc {
      template class RFieldDftComparison<1, CppTp<1> >;
      template class RFieldDftComparison<2, CppTp<2> >;
      template class RFieldDftComparison<3, CppTp<3> >;
   }
}
