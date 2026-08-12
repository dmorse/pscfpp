/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFields.h"             // class header
#include <cp/field/CFields.tpp>  // base class template implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Cp {
      template class CFields<1, Prdc::CField<1,CPT>, Cpc::FieldIo<1> >;
      template class CFields<2, Prdc::CField<2,CPT>, Cpc::FieldIo<2> >;
      template class CFields<3, Prdc::CField<3,CPT>, Cpc::FieldIo<3> >;
   }
   namespace Cpc {
      template class CFields<1>;
      template class CFields<2>;
      template class CFields<3>;
   }
}
