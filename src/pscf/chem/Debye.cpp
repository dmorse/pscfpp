/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Debye.h"
#include <cmath>

namespace Pscf {
namespace Debye {

   /*
   * Intrablock correlation function (thread model)
   */
   double dt(double ksq, double length, double kuhn)
   {
      double x = length*ksq*kuhn*kuhn/6.0;
      double d;
      if (x < 1.0E-5) {
         d = 1.0 - ((x*x)/3.0);
      } else {
         d = 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
      d *= length*length;
      return d;
   }

   /*
   * Intrablock correlation function (bead model)
   */
   double db(double ksq, double nBead, double kuhn)
   {
      double y = ksq*kuhn*kuhn/6.0;
      double x = nBead*y;
      double d;
      if (x < 1.0E-5) {
         d = 1.0 - (x*x)/3.0;
         d *= nBead*nBead;
      } else {
         double z = std::exp(-y);
         double u = 1 - z;
         d = 2.0 * z * (std::exp(-x) + nBead*u - 1.0) / (u*u);
         d += nBead;
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
   double eb(double ksq, double nBead, double kuhn)
   {
      double y = ksq*kuhn*kuhn/6.0;
      double x = nBead*y;
      double e;
      if (x < 1.0E-5) {
         e = 1.0 - x/2.0 + x*x/6.0;
         e *= nBead;
      } else {
         double z = std::exp(-y);
         e = (1.0 - std::exp(-x)) / (1 - z);
      }
      return e;
   }

} // namespace Debye
} // namespace Pscf
