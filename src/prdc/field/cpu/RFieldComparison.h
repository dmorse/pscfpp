#ifndef PRDC_CPU_R_FIELD_COMPARISON_H
#define PRDC_CPU_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/CppTp.h>             // specialized argument
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
   class RFieldComparison<D, CppTp<D> >
    : public FieldComparison< RField<D, CppTp<D> > >
   {};

   // Explicit instantiation declarations
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;

} // namespace Prdc
} // namespace Pscf
#endif
