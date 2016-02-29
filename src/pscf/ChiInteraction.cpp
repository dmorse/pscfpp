/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiInteraction.h"

namespace Pscf {
   
   using namespace Util;

   ChiInteraction::ChiInteraction()
    : Interaction()
   {  setClassName("ChiInteraction"); }

   ChiInteraction::~ChiInteraction()
   {}

   void ChiInteraction::readParameters(std::istream& in)
   {
      UTIL_CHECK(nMonomer() > 0);
      chi_.allocate(nMonomer(), nMonomer());
      readDSymmMatrix(in, "chi", chi_, nMonomer());
   }

   void 
   ChiInteraction::computeW(Array<double> const & c,
                                double p, Array<double>& w)
   {
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         w[i] = p;
         for (j = 0; j < nMonomer(); ++j) {
            w[i] += chi_(i, j)* c[j];
         }
      }
   }

   void 
   ChiInteraction::computeC(Array<double> const & w, Array<double>& c)
   {}

} // namespace Pscf
