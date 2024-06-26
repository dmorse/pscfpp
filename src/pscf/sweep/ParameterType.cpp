/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParameterType.h"

namespace Pscf {

   // Constructor
   ParameterType::ParameterType()
    : name(),
      nId(0),
      modifierPtr_(0)
   {}

   // Destructor
   ParameterType::~ParameterType()
   {}

} // namespace Pscf