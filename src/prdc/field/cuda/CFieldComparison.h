#ifndef PRDC_C_FIELD_COMPARISON_CU_H
#define PRDC_C_FIELD_COMPARISON_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/CUT.h>         // specialized template argument
#include <prdc/field/cuda/CField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class CFieldComparison;

   /**
   * Comparator for CField (k-grid) arrays, defined in device memory.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class CFieldComparison<D,CUT>
   {

   public:

      /**
      * Default constructor.
      *
      * Initializes maxDiff and rmsDiff to zero.
      */
      CFieldComparison();

      // Use compiler defined destructor and assignment operator.

      /**
      * Compare individual fields.
      *
      * Array dimensions must agree. An Exception is thrown if the 
      * capacities of fields a and b are not equal.
      *
      * \param a  1st field
      * \param b  2nd field
      * \return   maximum element-by-element difference (maxDiff)
      */ 
      double compare(CField<D,CUT> const & a, CField<D,CUT> const & b);

      /**
      * Compare arrays of fields associated with different monomer types.
      *
      * All array dimensions must agree.
      *
      * An exception is thrown if the capacities of the enclosing 
      * DArrays (the number of monomers) are not equal or if the
      * capacities of any pair of individual fields  (number of grid 
      * points or basis functions) are not equal.
      *
      * \param a  1st DArray of fields
      * \param b  2nd DArray of fields
      * \return   maximum element-by-element difference (maxDiff)
      */ 
      double compare(DArray< CField<D,CUT> > const & a, 
                     DArray< CField<D,CUT> > const & b);

      /**
      * Return the precomputed maximum element-by-element difference.
      *
      * This function returns the maximum difference between corresponding
      * field array elements found by the most recent comparison.
      */
      double maxDiff() const
      {  return maxDiff_; }
   
      /**
      * Return the precomputed root-mean-squared difference.
      *
      * This function returns the root-mean-squared difference between 
      * corresponding elements found by the most recent comparison.
      */
      double rmsDiff() const
      {  return rmsDiff_; }
   
   private:
  
      // Maximum element-by-element difference. 
      double maxDiff_;

      // Room-mean-squared element-by-element difference. 
      double rmsDiff_;
   
   };

   // Explicit instantiation declarations
   extern template class CFieldComparison<1,CUT>;
   extern template class CFieldComparison<2,CUT>;
   extern template class CFieldComparison<3,CUT>;


} // namespace Prdc
} // namespace Pscf
#endif
