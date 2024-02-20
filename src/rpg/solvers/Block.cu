/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
