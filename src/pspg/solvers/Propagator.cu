/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.tpp"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class Propagator<1>;
   template class Propagator<2>;
   template class Propagator<3>;

}
}
