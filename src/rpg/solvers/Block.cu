/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.tpp"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class Block<1>;
   template class Block<2>;
   template class Block<3>;

}
}
