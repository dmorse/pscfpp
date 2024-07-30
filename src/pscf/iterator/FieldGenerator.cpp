/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldGenerator.h"

namespace Pscf
{
   
   // Constructor
   FieldGenerator::FieldGenerator()
    : type_(None)
   {}

   // Destructor
   FieldGenerator::~FieldGenerator()
   {}

   // Allocate, check compatibility, calculate, and store the field(s)
   void FieldGenerator::setup()
   {
      allocate();
      checkCompatibility();
      if (!isGenerated()) {
         generate();
      } else {
         update();
      }
   }

   // Check whether system has changed and update the field(s) if necessary
   void FieldGenerator::update()
   {
      bool needed = updateNeeded();
      if (!needed) {
         // update not needed, do nothing
         return;
      } else {
         generate();
      }
   }   
   
}