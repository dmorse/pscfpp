/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixAndMatchEnvs.h"

namespace Pscf {
namespace Rpg {

   // FilmEnvironment instantiation
   template class FilmEnvironment<1>;
   template class FilmEnvironment<2>;
   template class FilmEnvironment<3>;

}
}