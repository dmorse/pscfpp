/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/SpaceGroup.tpp>

namespace Pscf {
namespace Prdc {

   template class SpaceGroup<1>;
   template class SpaceGroup<2>;
   template class SpaceGroup<3>;

   template void readGroup(std::string, SpaceGroup<1>& );
   template void readGroup(std::string, SpaceGroup<2>& );
   template void readGroup(std::string, SpaceGroup<3>& );

   template void writeGroup(std::string, SpaceGroup<1> const&);
   template void writeGroup(std::string, SpaceGroup<2> const&);
   template void writeGroup(std::string, SpaceGroup<3> const&);

}
}
