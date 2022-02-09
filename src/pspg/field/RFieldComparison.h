#ifndef PSPG_R_FIELD_COMPARISON_H
#define PSPG_R_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>
#include "RDField.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Comparator for fields in real-space (r-grid) format.
   * 
   * \ingroup Pspg_Field_Module
   */
   template <int D>
   class RFieldComparison
   {
   public:
      // Constructor
      RFieldComparison();

      // Comparator for individual fields. Note:challenging,because
      double compare(RDField<D> const& a, RDField<D> const& b);

      // Comparator for array of fields
      double compare(DArray<RDField<D>> const& a, DArray<RDField<D>> const& b);

      // Get maxDiff from FieldComparison
      double maxDiff() const
      {  return fieldComparison_.maxDiff(); }
      
      // Get rmsDiff from FieldComparison
      double rmsDiff() const
      {  return fieldComparison_.rmsDiff(); }
      

   private:

      // True if a comparison has been made, false otherwise.
      bool compared_;

      // Composition usage of FieldComparison, rather than inheritance.
      FieldComparison< DArray< cudaReal > > fieldComparison_;

   };


   #ifndef PSPG_R_FIELD_COMPARISON_TPP
   extern template class RFieldComparison<1>;
   extern template class RFieldComparison<2>;
   extern template class RFieldComparison<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
#endif
