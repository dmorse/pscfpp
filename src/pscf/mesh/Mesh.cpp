/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
