#ifndef PRDC_CUDA_C_FIELD_COMPARISON_H
#define PRDC_CUDA_C_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/CField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Comparator for CField (k-grid) arrays, defined in device memory.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class CFieldComparison {

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
      double compare(CField<D> const & a, CField<D> const & b);

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
      double compare(DArray< CField<D> > const & a, 
                     DArray< CField<D> > const & b);

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

   #ifndef PRDC_CUDA_C_FIELD_COMPARISON_TPP
   // Suppress implicit instantiation
   extern template class CFieldComparison<1>;
   extern template class CFieldComparison<2>;
   extern template class CFieldComparison<3>;
   #endif


} // namespace Prdc::Cuda
} // namespace Prdc
} // namespace Pscf
#endif
