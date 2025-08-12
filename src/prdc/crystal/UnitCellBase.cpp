/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellBase.tpp"

#include <util/math/Constants.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   template class UnitCellBase<1>;
   template class UnitCellBase<2>;
   template class UnitCellBase<3>;

}
}
