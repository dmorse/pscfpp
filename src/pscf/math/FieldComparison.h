#ifndef PSCF_FIELD_COMPARISON_H
#define PSCF_FIELD_COMPARISON_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace Pscf {

   using namespace Util;

   /**
   * Comparison of element-by-element differences between field arrays.
   *
   * The template argument FT may be RField<D> for representations of a
   * field or fields on an r-grid or DArray<double> for representation 
   * using a symmetry-adapted basis.
   */
   template <class FT>
   class FieldComparison {

   public:

      /**
      * Default constructor.
      *
      * Initializes maxDiff and rmsDiff to zero.
      */
      FieldComparison(int begin = 0);

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
      double compare(FT const& a, FT const& b);

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
      double compare(DArray<FT> const& a, DArray<FT> const& b);

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
   
      // Index of first element (0 or 1)
      int begin_;
   
   };

} // namespace Pscf
#include "FieldComparison.tpp"
#endif
