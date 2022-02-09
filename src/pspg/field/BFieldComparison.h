#ifndef PSPC_B_FIELD_COMPARISON_H
#define PSPC_B_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>
//#include <pspg/field/RDField.h>
#include <util/containers/DArray.h>
#include <pspg/GpuResources.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Comparator for fields in symmetry-adapted basis format.
   */
   class BFieldComparison : public FieldComparison< DArray <cudaReal*> >
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

      double compare(DArray<cudaReal*> const& a, DArray<cudaReal*> const& b, int nStar);

   private:

   };

} // namespace Pspc
} // namespace Pscf
#endif
