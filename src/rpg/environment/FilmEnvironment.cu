/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/environment/FilmEnvironment.h>

#include <rpg/environment/FilmFieldGenMask.h>
#include <rpg/environment/FilmFieldGenExt.h>

#include <rp/environment/FilmEnvironment.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FilmEnvironment<1, Rpg::Types<1> >;
      template class FilmEnvironment<2, Rpg::Types<2> >;
      template class FilmEnvironment<3, Rpg::Types<3> >;
   }
}
