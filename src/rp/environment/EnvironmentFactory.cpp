/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/environment/EnvironmentFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EnvironmentFactory<1,CPT>;
      template class EnvironmentFactory<2,CPT>;
      template class EnvironmentFactory<3,CPT>;
   }
}
