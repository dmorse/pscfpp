#ifndef PSCF_B_FIELD_COMPARISON_H
#define PSCF_B_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/FieldComparison.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Comparator for fields in symmetry-adapted basis format.
   * 
   * \ingroup Pscf_Prdc_Crystal_Module
   */
   class BFieldComparison : public FieldComparison< DArray<double> >
   {

   public:

      /**
      * Constructor.
      * 
      * The basis function with index 0 in a symmetry adapted basis is 
      * always a spatially homogeneous function, i.e., a constant. In 
      * some situations, we may be interested in determining whether
      * two fields are equivalent to within a constant.
      *
      * Set begin = 0, which is the default, to include the coefficient 
      * of the first basis function in the comparison, thus determining
      * how close to fields are to being strictly equal. 
      *
      * Set begin = 1 to exclude the coefficient of the first basis 
      * function, thus comparing only deviatoric parts of the fields.
      *
      * \param begin  index of first element to include in comparison.
      */
      BFieldComparison(int begin = 0);

   };

} // namespace Prdc
} // namespace Pscf
#endif
