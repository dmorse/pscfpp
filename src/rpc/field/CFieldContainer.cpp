/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CFieldContainer.h"
#include <rpc/field/FieldIo.h>
#include <prdc/cpu/RField.h>
#include <prdc/field/CFieldsTmpl.tpp>

namespace Pscf {
   namespace Prdc {
      // Explicit instantiations of base class template
      template class CFieldsTmpl<1, Cpu::RField<1>, Rpc::FieldIo<1> >;
      template class CFieldsTmpl<2, Cpu::RField<2>, Rpc::FieldIo<2> >;
      template class CFieldsTmpl<3, Cpu::RField<3>, Rpc::FieldIo<3> >;
   } 
   namespace Rpc {
      // Explicit instantiations of this class
      template class CFieldContainer<1>;
      template class CFieldContainer<2>;
      template class CFieldContainer<3>;
   
   } 
}
