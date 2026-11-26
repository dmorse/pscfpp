/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldGeneratorBase.h"

namespace Pscf
{
   
   // Constructor
   FieldGeneratorBase::FieldGeneratorBase()
    : type_(None),
      isDependent_(false)
   {}

   // Destructor
   FieldGeneratorBase::~FieldGeneratorBase()
   {}

   // Checks compatibility, and (re)generates the fields if needed.
   void FieldGeneratorBase::generate()
   {
      // Check the system even if the fields need no update: a change
      // that leaves the fields unchanged can make it incompatible
      checkCompatibility();

      if (needsUpdate()) {
         compute();
      }
   }   
   
}