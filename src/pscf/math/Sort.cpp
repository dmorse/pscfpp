/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sort.h"
#include <util/global.h>
#include <vector>
#include <algorithm>

namespace Pscf {
namespace Sort {

   using namespace Util;

   /*
   * Sort a std::vector of items.
   */
   template <typename T>
   void sort(std::vector< Item<T> >& items)
   {  std::sort(items.begin(), items.end()); }

   /*
   * Identify slices of equal values within a sorted vector.
   */
   template <typename T>
   void findBunches(std::vector< Item<T> > const & items,
                    std::vector< Slice >& bunches,
                    T epsilon)
   {
      int size = items.size();
      UTIL_CHECK(size > 0);
      bunches.clear();
      Slice bunch;
      bunch.begin = 0;
      bunch.end = 0;
      T newVal = items[0].value;
      T oldVal = newVal;
      if (size > 1) {
         for (int i = 1; i < size; ++i) {
            newVal = items[i].value;
            UTIL_CHECK(newVal > oldVal - epsilon);
            if (newVal > oldVal + epsilon) {
               bunch.end = i;
               UTIL_CHECK(bunch.end > bunch.begin);
               bunches.push_back(bunch);
               bunch.begin = i;
               bunch.end = i;
               oldVal = newVal;
            }
         }
      }
      UTIL_CHECK(bunch.begin < size);
      bunch.end = size;
      UTIL_CHECK(bunch.end > bunch.begin);
      bunches.push_back(bunch);
   }

   // Explicit instantiation definitions

   template void sort<double>(std::vector< Item<double> >& );
   template void sort<float>(std::vector< Item<float> >& );

   template
   void findBunches<double>(std::vector< Item<double> > const &,
                            std::vector< Slice >&, double);

   template
   void findBunches<float>(std::vector< Item<float> > const &,
                           std::vector< Slice >&, float);

} // namespace Sort
} // namespace Pscf
