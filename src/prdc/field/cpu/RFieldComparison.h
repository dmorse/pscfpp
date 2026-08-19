#ifndef PRDC_CPU_R_FIELD_COMPARISON_H
#define PRDC_CPU_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>             // specialized argument
#include <pscf/math/FieldComparison.h>  // base class template
#include "RField.h"                     // base class template argument

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class RFieldComparison;

   /**
   * Comparator for fields in real-space (r-grid) format.
   * 
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class RFieldComparison<D,CPT>
    : public FieldComparison< RField<D,CPT> >
   {};

   // Explicit instantiation declarations
   extern template class RFieldComparison<1,CPT>;
   extern template class RFieldComparison<2,CPT>;
   extern template class RFieldComparison<3,CPT>;

} // namespace Prdc
} // namespace Pscf
#endif
