/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParameterType.h"

namespace Pscf {

   // Constructor
   ParameterType::ParameterType()
    : name(),
      nId(0),
      modifierPtr(0)
   {}

   // Alternate constructor that sets all members
   ParameterType::ParameterType(std::string name, int nId, 
                                    ParameterModifier& modifier)
    : name(name),
      nId(nId),
      modifierPtr(&modifier)
   {}

   // Destructor
   ParameterType::~ParameterType()
   {}

} // namespace Pscf