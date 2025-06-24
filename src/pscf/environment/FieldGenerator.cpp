/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldGenerator.h"

namespace Pscf
{
   
   // Constructor
   FieldGenerator::FieldGenerator()
    : type_(None),
      isDependent_(false)
   {}

   // Destructor
   FieldGenerator::~FieldGenerator()
   {}

   // Checks if fields need to be (re)generated. If so, generates them. 
   void FieldGenerator::generate()
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