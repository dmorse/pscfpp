/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDftComparison.tpp"

namespace Pscf {
   namespace Prdc {
      template class RFieldDftComparison<1,CPT>;
      template class RFieldDftComparison<2,CPT>;
      template class RFieldDftComparison<3,CPT>;
   }
}
