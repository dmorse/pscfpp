/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/FilmFieldGenMask_u.h>
#include <rp/environment/FilmFieldGenExt_u.h>

#include <pscf/backend/CUT.h>
#include <rp/environment/FilmEnvironment.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FilmEnvironment<1,CUT>;
      template class FilmEnvironment<2,CUT>;
      template class FilmEnvironment<3,CUT>;
   }
}
