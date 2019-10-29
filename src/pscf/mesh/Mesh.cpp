/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mesh.h"
#include <util/global.h>

namespace Pscf
{

   using namespace Util;

   template class Mesh<1>;
   template class Mesh<2>;
   template class Mesh<3>;

}
