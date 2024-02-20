#ifndef PSPC_B_FIELD_COMPARISON_H
#define PSPC_B_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>
//#include <rpg/field/RDField.h>
#include <util/containers/DArray.h>
#include <rpg/field/DField.h>
#include <rpg/math/GpuResources.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Comparator for fields in symmetry-adapted basis format.
   */
   class BFieldComparison
   {

   public:

      /**
      * Constructor.
      * 
      * The basis function with index 0 in a symmetry adapted basis is always
      * a spatially homogeneous function, i.e., a constant. In some situations,
      * we may be interested in determining whether two fields are equivalent
      * to within a constant.
      *
      * Set begin = 0 to include the coefficient of the first basis function
      * in the comparison, thus determining how close to fields are to being
      * strictly equal. 
      *
      * Set begin = 1 to exclude the coefficient of the first basis function, 
      * thus comparing only deviatoric parts of the fields.
      *
      * \param begin  index of first element to include in comparison.
      */
      BFieldComparison(int begin = 0);

      double compare(DField<cudaReal> const& a, DField<cudaReal> const& b);

      double compare(DArray<DField<cudaReal>> const& a, DArray<DField<cudaReal>> const& b);

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

} // namespace Pspc
} // namespace Pscf
#endif
