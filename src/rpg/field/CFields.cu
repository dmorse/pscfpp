/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"              // class header
#include <rp/field/CFields.tpp>   // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      template class CFields<1, RField<1>, Rpg::FieldIo<1> >;
      template class CFields<2, RField<2>, Rpg::FieldIo<2> >;
      template class CFields<3, RField<3>, Rpg::FieldIo<3> >;
   } 
   namespace Rpg {
      template class CFields<1>;
      template class CFields<2>;
      template class CFields<3>;
   } 
}
