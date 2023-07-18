/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.tpp"

namespace Pscf {
namespace Pspc {

   template<> long Analyzer<1>::baseInterval = 0;
   template<> long Analyzer<2>::baseInterval = 0;
   template<> long Analyzer<3>::baseInterval = 0;

   template class Analyzer<1>;
   template class Analyzer<2>;
   template class Analyzer<3>;
}
}
