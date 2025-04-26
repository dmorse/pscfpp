/*
* PSCF Package
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldBasisConverter.tpp"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   // Explicit class instantiations
   template class FieldBasisConverter<1>;
   template class FieldBasisConverter<2>;
   template class FieldBasisConverter<3>;

}
}
}
