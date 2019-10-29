/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
