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
namespace Rp {

   // Explicit instantiation declarations
   template class FilmFieldGenMask<1, Rpg::Types<1> >;
   template class FilmFieldGenMask<2, Rpg::Types<2> >;
   template class FilmFieldGenMask<3, Rpg::Types<3> >;
   
}
}
