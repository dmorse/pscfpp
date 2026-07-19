/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/environment/FilmEnvironment.h>

#include <rpc/environment/FilmFieldGenMask.h>
#include <rpc/environment/FilmFieldGenExt.h>

#include <rp/environment/FilmEnvironment.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FilmEnvironment<1, Cpp<1> >;
      template class FilmEnvironment<2, Cpp<2> >;
      template class FilmEnvironment<3, Cpp<3> >;
   }
}
