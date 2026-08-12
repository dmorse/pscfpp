/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/PerturbationDerivative.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationDerivative<1,CPT>;
      template class PerturbationDerivative<2,CPT>;
      template class PerturbationDerivative<3,CPT>;
   }
}
