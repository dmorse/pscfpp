#ifndef PSPC_R_FIELD_COMPARISON_H
#define PSPC_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/field/FieldComparison.h>
#include "RField.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Comparator for fields in real-space (r-grid) format.
   */
   template <int D>
   class RFieldComparison : public FieldComparison< RField<D> >
   {};

   #ifndef PSPC_R_FIELD_COMPARISON_CPP
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
