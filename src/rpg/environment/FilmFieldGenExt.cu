/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenExt.tpp"
#include <rpg/field/Domain.cu>
#include <rpg/field/CFields.cu>
#include <rpg/field/WFields.cu>
#include <rpg/field/Mask.cu>

namespace Pscf {
namespace Rpg {

   template class FilmFieldGenExt<1>;
   template class FilmFieldGenExt<2>;
   template class FilmFieldGenExt<3>;

}
}
