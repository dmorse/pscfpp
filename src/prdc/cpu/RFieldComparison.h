#ifndef PRDC_CUDA_R_FIELD_COMPARISON_H
#define PRDC_CUDA_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>
#include "RField.h"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Comparator for fields in real-space (r-grid) format.
   * 
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class RFieldComparison : public FieldComparison< RField<D> >
   {};

   #ifndef PRDC_R_FIELD_COMPARISON_CPP
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;
   #endif

} // namespace Pscf::Prdc::Gpu
} // namespace Pscf::Prdc
} // namespace Pscf
#endif
