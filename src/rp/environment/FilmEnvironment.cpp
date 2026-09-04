/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/FilmFieldGenMask_c.h>
#include <rp/environment/FilmFieldGenExt_c.h>

#include <pscf/backend/cpp/CPT.h>
#include <rp/environment/FilmEnvironment.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FilmEnvironment<1,CPT>;
      template class FilmEnvironment<2,CPT>;
      template class FilmEnvironment<3,CPT>;
   }
}
