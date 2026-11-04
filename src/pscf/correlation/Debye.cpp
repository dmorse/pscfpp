/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Debye.h"
#include <cmath>

namespace Pscf {
namespace Correlation {

   /*
   * Intrablock correlation function (thread model)
   */
   double dt(double ksq, double length, double kuhn)
   {
      double x = length*ksq*kuhn*kuhn/6.0;
      double d;
      if (x < 1.0E-5) {
         d = 1.0 - (x/3.0) + (x*x/12.0);
      } else {
         d = 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
      d *= length*length;
      return d;
   }

   /*
   * Intrablock correlation function (bead model)
   */
   double db(double ksq, double n, double kuhn)
   {
      double y = ksq*kuhn*kuhn/6.0;
      double x = n*y;
      double d;
      if (x < 1.0E-5) {
         double nSq = n*n;
         d = nSq - x*(nSq - 1.0)/3.0 + x*x*(nSq - 1.0)/12.0;
      } else {
         const double z = std::exp(-y);
         const double u = 1 - z;
         d = 2.0 * z * (std::exp(-x) + n*u - 1.0) / (u*u);
         d += n;
      }
      return d;
   }

   /*
   * End factor for one block (thread model)
   */
   double et(double ksq, double length, double kuhn)
   {
      double x = ksq*length*kuhn*kuhn/6.0;
      double e;
      if (x < 1.0E-5) {
         e = 1.0 - x/2.0 + x*x/6.0;
      } else {
         e = (1.0 - std::exp(-x)) / x;
      }
      e *= length;
      return e;
   }

   /*
   * End factor for one block (bead model)
   */
   double eb(double ksq, double n, double kuhn)
   {
      double y = ksq*kuhn*kuhn/6.0;
      double x = n*y;
      double e;
      if (x < 1.0E-5) {
         e = n - x*(n - 1.0)/2.0;
         e += x*y*(n*n - 1.5*n + 0.5)/6.0;
      } else {
         const double z = std::exp(-y);
         e = (1.0 - std::exp(-n*y)) / (1 - z);
      }
      return e;
   }

} // namespace Correlation
} // namespace Pscf
