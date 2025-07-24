/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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

   // Checks if fields need to be (re)generated. If so, generates them. 
   void FieldGeneratorBase::generate()
   {
      if (needsUpdate()) {
         checkCompatibility();
         compute();
      } else {
         // update not needed, do nothing
         return;
      }
   }   
   
}