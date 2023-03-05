/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorOld.tpp"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class AmIteratorOld<1>;
   template class AmIteratorOld<2>;
   template class AmIteratorOld<3>;

}
}
