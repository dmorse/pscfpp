#ifndef PSCF_WAVE_NORM_H
#define PSCF_WAVE_NORM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/Sort.h>
#include <util/containers/GArray.h>
#include <util/global.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace Pscf {

   /**
   * Are two squared wavevector norms equal?
   *
   * Two values are taken to be equal if they differ by at most 1.0E-8,
   * or by at most 1.0E-10 relative to the larger one. The relative
   * tolerance allows for rounding errors of large squared norms (e.g.,
   * in small unit cells), which can exceed 1.0E-8.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param a  squared norm of a wavevector
   * \param b  squared norm of another wavevector
   */
   inline bool equalWaveNorm(double a, double b)
   {
      double diff = std::abs(a - b);
      return (diff <= 1.0E-8) || (diff <= 1.0E-10*std::max(a, b));
   }

   /**
   * Identify bunches of equal squared wavevector norms.
   *
   * Given items sorted by increasing squared norm (Sort::sort), find
   * the slices of items whose values are equal to that of the first
   * item of the slice, as defined by equalWaveNorm.
   *
   * \ingroup Pscf_Crystal_Module
   *
   * \param items  items sorted by increasing squared norm (in)
   * \param bunches  slices of items of equal squared norm (out)
   */
   inline
   void findWaveBunches(std::vector< Sort::Item<double> > const & items,
                        GArray< Sort::Bunch >& bunches)
   {
      int size = items.size();
      UTIL_CHECK(size > 0);
      bunches.clear();
      Sort::Bunch bunch;
      bunch[0] = 0;
      for (int i = 1; i < size; ++i) {
         if (!equalWaveNorm(items[i].value, items[bunch[0]].value)) {
            bunch[1] = i;
            bunches.append(bunch);
            bunch[0] = i;
         }
      }
      bunch[1] = size;
      bunches.append(bunch);
   }

}
#endif
