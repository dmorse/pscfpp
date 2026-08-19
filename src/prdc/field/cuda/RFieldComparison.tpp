#ifndef PRDC_R_FIELD_COMPARISON_CU_TPP
#define PRDC_R_FIELD_COMPARISON_CU_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.h"

#include <pscf/cuda/HostDArray.h>

namespace Pscf {
namespace Prdc {

   /*
   * Default constructor.
   */
   template <int D>
   RFieldComparison<D,CUT>::RFieldComparison()
   {};

   /*
   * Comparator for individual fields.
   */
   template <int D>
   double RFieldComparison<D,CUT>::compare(RField<D,CUT> const& a, 
                                       RField<D,CUT> const& b)
   {
      int nPoints = a.capacity();

      // Copy fields a,b to local arrays ha, hb on the host CPU
      HostDArray<cudaReal> ha(nPoints);
      HostDArray<cudaReal> hb(nPoints);
      ha = a;
      hb = b;

      fieldComparison_.compare(ha, hb);
      compared_ = true;

      return fieldComparison_.maxDiff();
   }

   /*
   * Comparator for arrays of fields.
   */
   template <int D>
   double RFieldComparison<D,CUT>::compare(DArray< RField<D,CUT> > const& a, 
                                       DArray< RField<D,CUT> > const& b)
   {
      int nFields = a.capacity();
      int nPoints = a[0].capacity();

      // Copy fields to HostDArray containers on CPU host
      DArray< HostDArray<cudaReal> > ha, hb;
      ha.allocate(nFields);
      hb.allocate(nFields);
      for (int i = 0; i < nFields; i++) {
         ha[i] = a[i];
         hb[i] = b[i];
      }

      // Perform comparison
      fieldComparison_.compare(ha, hb);
      compared_ = true;

      return fieldComparison_.maxDiff();
   }

} // namespace Prdc
} // namespace Pscf
#endif
