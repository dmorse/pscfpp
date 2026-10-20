#ifndef PSCF_CPP_R_FIELD_COMPARISON_H
#define PSCF_CPP_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>  // base class template
#include <pscf/backend/cpp/CPT.h>       // backend identifier class
#include <pscf/field/cpp/RField.h>      // template argument

namespace Pscf {

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

} // namespace Pscf
#endif
