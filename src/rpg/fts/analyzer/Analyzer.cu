/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.tpp"

namespace Pscf {
namespace Rpg {

   template<> long Analyzer<1>::baseInterval = 1;
   template<> long Analyzer<2>::baseInterval = 1;
   template<> long Analyzer<3>::baseInterval = 1;

   template class Analyzer<1>;
   template class Analyzer<2>;
   template class Analyzer<3>;
}
}
