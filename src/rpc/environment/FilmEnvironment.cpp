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
      template class FilmEnvironment<1,CPT>;
      template class FilmEnvironment<2,CPT>;
      template class FilmEnvironment<3,CPT>;
   }
}
