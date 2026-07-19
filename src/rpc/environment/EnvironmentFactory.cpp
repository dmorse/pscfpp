/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/environment/EnvironmentFactory.h>  // header

#include "FilmEnvironment.h"

#include <rp/environment/EnvironmentFactory.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EnvironmentFactory<1, Rpc::Types<1> >;
      template class EnvironmentFactory<2, Rpc::Types<2> >;
      template class EnvironmentFactory<3, Rpc::Types<3> >;
   }
}
