#ifndef PRDC_CUDA_R_FIELD_COMPARISON_H
#define PRDC_CUDA_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include <pscf/math/FieldComparison.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Comparator for fields in real-space (r-grid) format.
   * 
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class RFieldComparison
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
      double compare(RField<D> const& a, RField<D> const& b);

      /**
      * Comparator for arrays of fields.
      *
      * \param a first array of fields
      * \param b second array of fields
      */
      double 
      compare(DArray< RField<D> > const& a, DArray< RField<D> > const& b);

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
      FieldComparison< HostDArray< cudaReal > > fieldComparison_;

   };

   #ifndef PRDC_CUDA_R_FIELD_COMPARISON_TPP
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;
   #endif

} // namespace Prdc::Cuda
} // namespace Prdc
} // namespace Pscf
#endif
