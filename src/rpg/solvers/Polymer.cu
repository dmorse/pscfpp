/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.tpp"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class Polymer<1>;
   template class Polymer<2>;
   template class Polymer<3>;

}
}
