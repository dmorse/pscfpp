/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SymmetryGroup.tpp>
#include <pscf/crystal/SpaceSymmetry.tpp>

namespace Pscf
{

   template class SymmetryGroup< SpaceSymmetry<1> >;
   template class SymmetryGroup< SpaceSymmetry<2> >;
   template class SymmetryGroup< SpaceSymmetry<3> >;

}
