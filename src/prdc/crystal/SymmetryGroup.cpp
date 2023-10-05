/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/SymmetryGroup.tpp>
#include <prdc/crystal/SpaceSymmetry.tpp>

namespace Pscf {
namespace Prdc {

   template class SymmetryGroup< SpaceSymmetry<1> >;
   template class SymmetryGroup< SpaceSymmetry<2> >;
   template class SymmetryGroup< SpaceSymmetry<3> >;

}
}
