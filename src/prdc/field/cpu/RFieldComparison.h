#ifndef PRDC_CPU_R_FIELD_COMPARISON_H
#define PRDC_CPU_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>  // base class template
#include "RField.h"                     // base class template argument

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /**
   * Comparator for fields in real-space (r-grid) format.
   * 
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class RFieldComparison : public FieldComparison< RField<D> >
   {};

   // Explicit instantiation declarations
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;

} // namespace Cpu
} // namespace Prdc
} // namespace Pscf
#endif
