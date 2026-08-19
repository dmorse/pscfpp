/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenExtBase.tpp"

// Explicit instantiation definitions
namespace Pscf {
   namespace Prdc {
      template class FilmFieldGenExtBase<1>;
      template class FilmFieldGenExtBase<2>;
      template class FilmFieldGenExtBase<3>;
   }
}
