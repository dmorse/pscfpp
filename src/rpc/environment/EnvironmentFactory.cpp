/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/EnvironmentFactory.h>  // header

#include <rp/environment/FilmEnvironment.h>

#include <rp/environment/EnvironmentFactory.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class EnvironmentFactory<1,CPT>;
      template class EnvironmentFactory<2,CPT>;
      template class EnvironmentFactory<3,CPT>;
   }
}
