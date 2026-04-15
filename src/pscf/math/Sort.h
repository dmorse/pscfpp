#ifndef PSCF_SORT_H
#define PSCF_SORT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <vector>

namespace Pscf {

   /**
   * Utilities for sorting real values.
   *
   * \defgroup Pscf_Math_Sort_Module
   * \ingroup Pscf_Math_Module
   */
   namespace Sort {

      /**
      * Struct with value and index, to keep track of permutation.
      *
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      struct Item {
         T value;
         int id;
      };

      /**
      * Bounds of a contiguous array slice.
      *
      * \ingroup Pscf_Math_Sort_Module
      */
      struct Slice {
         int begin;
         int end;
      };

      /**
      * Less than comparator for Item<T> objects.
      *
      * Comparison based on Item<T>::value member.
      *
      * \param a  1st item
      * \param b  2nd item
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      bool operator < (Item<T> const & a, Item<T> const & b)
      {  return (a.value < b.value); }

      /**
      * Sort a std::vector< Item<T> > by ascending item value.
      *
      * On return, array "items" is sorted with non-descending values.
      *
      * \param items  array of comparable items
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      void sort(std::vector< Item<T> >& items);

      /**
      * Identify "bunches" of equal values within a sorted vector.
      *
      * On entry, array "items" must be sorted on entry. The sort
      * function must thus normally be called before this.
      *
      * \param items  sorted array of items
      * \param bunches  array of slices of items with equal value
      * \param epsilon  tolerance
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      void findBunches(std::vector< Item<T> > const & items,
                       std::vector< Slice >& bunches,
                       T epsilon);


      // Explicit instantiation declarations

      extern template void sort<double>(std::vector< Item<double> >& );
      extern template void sort<float>(std::vector< Item<float> >& );

      extern template
      void findBunches<double>(std::vector< Item<double> > const &,
                               std::vector< Slice >&, double);

      extern template
      void findBunches<float>(std::vector< Item<float> > const &,
                              std::vector< Slice >&, float);

   } // namespace Sort

} // namespace Pscf
#endif
