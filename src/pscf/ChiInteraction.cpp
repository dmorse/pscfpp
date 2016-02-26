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
      readSymmDMatrix(in, "chi", chi_, nMonomer());
   }

} // namespace Pscf
