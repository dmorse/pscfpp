/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenMask.tpp"
#include <rpg/field/Domain.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>

namespace Pscf {
namespace Rpg {

   // Class declarations
   template class FilmFieldGenMask<1>;
   template class FilmFieldGenMask<2>;
   template class FilmFieldGenMask<3>;
   
}
}
