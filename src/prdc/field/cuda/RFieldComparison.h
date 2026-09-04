#ifndef PRDC_R_FIELD_COMPARISON_CU_H
#define PRDC_R_FIELD_COMPARISON_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/CUT.h>           // specialized argument
#include <pscf/math/FieldComparison.h>  // member
#include "RField.h"

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class RFieldComparison;

   /**
   * Comparator for fields in real-space (r-grid) format.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class RFieldComparison<D,CUT>
   {
   public:

      /**
      * Constructor.
      */
      RFieldComparison();

      /**
      * Comparator for individual fields.
      *
      * \param a first array of fields
      * \param b second array of fields
      */
      double compare(RField<D,CUT> const& a,
                     RField<D,CUT> const& b);

      /**
      * Comparator for arrays of fields.
      *
      * \param a first array of fields
      * \param b second array of fields
      */
      double
      compare(DArray< RField<D,CUT> > const& a,
              DArray< RField<D,CUT> > const& b);

      /**
      * Get precomputed maximum element-by-element difference.
      */
      double maxDiff() const
      {  return fieldComparison_.maxDiff(); }

      /**
      * Get precomputed rms difference.
      */
      double rmsDiff() const
      {  return fieldComparison_.rmsDiff(); }

   private:

      // True if a comparison has been made, false otherwise.
      bool compared_;

      // Use FieldComparison template via composition
      FieldComparison< HostDArray<cudaReal> > fieldComparison_;

   };

   extern template class RFieldComparison<1,CUT>;
   extern template class RFieldComparison<2,CUT>;
   extern template class RFieldComparison<3,CUT>;

} // namespace Prdc
} // namespace Pscf
#endif
