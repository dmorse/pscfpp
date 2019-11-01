/*
* PSCF Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.tpp"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class Mixture<1>;
   template class Mixture<2>;
   template class Mixture<3>;

}
}
