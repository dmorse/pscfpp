#ifndef UTIL_FEQ_H
#define UTIL_FEQ_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <cmath>

namespace Util
{

   /**
   * Are two floating point numbers equal to within round-off error?
   * 
   * \ingroup Math_Module
   */
   inline bool feq(double x, double y, double eps = 1.0E-10)
   {
      double diff = fabs(x - y) / (fabs(x) + fabs(y) + 1.0E-5);
      return (diff < eps);
   }

}
#endif
