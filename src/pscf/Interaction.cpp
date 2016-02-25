/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"

namespace Pscf {
   
   using namespace Util;

   Interaction::Interaction()
    : nMonomer_(0)
   {  setClassName("Interaction"); }

   Interaction::~Interaction()
   {}

   void Interaction::setNMonomer(int nMonomer)
   {  nMonomer_ = nMonomer; }

} // namespace Pscf
