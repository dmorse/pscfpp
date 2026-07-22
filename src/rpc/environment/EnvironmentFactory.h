#ifndef RPC_ENVIRONMENT_FACTORY_H
#define RPC_ENVIRONMENT_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/EnvironmentFactory.h>
#include <pscf/backends/CPT.h>

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class EnvironmentFactory<1,CPT>;
      extern template class EnvironmentFactory<2,CPT>;
      extern template class EnvironmentFactory<3,CPT>;
   }
}
#endif
