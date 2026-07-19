/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnvironmentFactory.h"

#include "FilmEnvironment.h"

#include <rp/environment/EnvironmentFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EnvironmentFactory<1, Rpg::Types<1> >;
      template class EnvironmentFactory<2, Rpg::Types<2> >;
      template class EnvironmentFactory<3, Rpg::Types<3> >;
   }
}
