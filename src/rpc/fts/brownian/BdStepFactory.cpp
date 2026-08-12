/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdStepFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BdStepFactory<1,CPT>;
      template class BdStepFactory<2,CPT>;
      template class BdStepFactory<3,CPT>;
   }
}
