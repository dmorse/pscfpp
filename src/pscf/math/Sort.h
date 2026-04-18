#ifndef PSCF_SORT_H
#define PSCF_SORT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Pair.h>
#include <util/containers/GArray.h>
#include <vector>

namespace Pscf {

   using namespace Util;

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
      * The sort function sorts Item objects in ascending order of
      * the "value" member variable. The id member variable is an
      * immutable integer identifier that is unchanged by sorting.
      *
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      struct Item {
         T value;
         int id;
      };

      /*
      * A Bunch represents a slice of the sorted array of items.
      *
      * A Bunch is used to represent a contiguous slice of the sorted
      * array of items produced by the Sort::sort function within which
      * values are equal to within some tolerance. If x is a Bunch, then
      * x[0] is the array index of the first element in a slice, and
      * x[1] an integer is one greater than the index of the last element
      * in that slice.
      *
      * \ingroup Pscf_Math_Sort_Module
      */
      using Bunch = Pair<int>;

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
      * The array "items" must be sorted on entry to this function. 
      * The "sort" function is thus called before this.
      *
      * On return, each element of array "bunches" contains the begin
      * and end indices delimiting a "bunch" of contiguous items in 
      * the sorted "items" array with equal value. Each item belongs 
      * to one bunch, and each bunch contains one or more items. 
      * Bunches are listed in order increasing beginning index, and 
      * span the entire array.
      *
      * \param items  sorted array of items
      * \param bunches  array of bunches of items with equal value
      * \param epsilon  tolerance
      * \ingroup Pscf_Math_Sort_Module
      */
      template <typename T>
      void findBunches(std::vector< Item<T> > const & items,
                       GArray< Bunch >& bunches,
                       T epsilon);


      // Explicit instantiation declarations

      extern template void sort<double>(std::vector< Item<double> >& );
      extern template void sort<float>(std::vector< Item<float> >& );

      extern template
      void findBunches<double>(std::vector< Item<double> > const &,
                               GArray< Bunch >&, double);

      extern template
      void findBunches<float>(std::vector< Item<float> > const &,
                              GArray< Bunch >&, float);

   } // namespace Sort

} // namespace Pscf
#endif
