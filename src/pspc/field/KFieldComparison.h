#ifndef PSPC_K_FIELD_COMPARISON_H
#define PSPC_K_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>
#include <pspc/field/RFieldDft.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /**
   * Comparator for RFieldDft (k-grid) arrays.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class KFieldComparison {

   public:

      /**
      * Default constructor.
      *
      * Initializes maxDiff and rmsDiff to zero.
      */
      KFieldComparison();

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
      double compare(RFieldDft<D> const& a, RFieldDft<D> const& b);

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
      * \param a  1st DArray of field
      * \param b  2nd DArray of field
      * \return   maximum element-by-element difference (maxDiff)
      */ 
      double compare(DArray<RFieldDft<D> > const& a, 
                     DArray<RFieldDft<D> > const& b);

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

   #ifndef PSPC_K_FIELD_COMPARISON_TPP
   // Suppress implicit instantiation
   extern template class KFieldComparison<1>;
   extern template class KFieldComparison<2>;
   extern template class KFieldComparison<3>;
   #endif


} // namespace Pspc
} // namespace Pscf
#endif
