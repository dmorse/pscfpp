/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mesh.tpp"
#include <util/global.h>

namespace Pscf
{

   using namespace Util;

   template class Mesh<1>;
   template class Mesh<2>;
   template class Mesh<3>;

}
