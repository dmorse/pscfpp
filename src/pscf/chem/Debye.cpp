/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Debye.h"
#include <cmath>

namespace Pscf {
namespace Debye {

   /*
   * Homopolymer correlation function.
   */
   double d(double ksq, double length, double kuhn)
   {
      double x = ksq*length*kuhn*kuhn/6.0;
      double g;
      if (x < 1.0E-5) {
         g = 1.0 - x*x/3.0;
      } else {
         g = 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
      return length*length*g;
   }

   /*
   * End factor for one block. 
   */
   double e(double ksq, double length, double kuhn)
   {
      double x = ksq*length*kuhn*kuhn/6.0;
      double h;
      if (x < 1.0E-5) {
         h = 1.0 - x/2.0 + x*x/6.0;
      } else {
         h = (1.0 - std::exp(-x)) / x;
      }
      return length*h;
   }

} // namespace Debye
} // namespace Pscf
