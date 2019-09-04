/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#define PSCF_SPACE_GROUP_CPP
#include <pscf/crystal/SpaceGroup.h>

namespace Pscf
{

   template class SpaceGroup<1>;
   template class SpaceGroup<2>;
   template class SpaceGroup<3>;

}
