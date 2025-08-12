/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldGenerator.h"

namespace Pscf {
namespace Prdc {

   // Constructor
   FieldGenerator::FieldGenerator()
   {}

   // Destructor
   FieldGenerator::~FieldGenerator()
   {}

   // Get contribution to the stress from this imposed field.
   double FieldGenerator::stress(int paramId) const
   // (default implementation)
   {  UTIL_THROW("Unimplemented stress() method called."); } 

   // Modify stress to minimize a property other than fHelmholtz. 
   // (default implementation)
   double FieldGenerator::modifyStress(int paramId, double stress) 
   const
   {  return stress; }

} // namespace Prdc
} // namespace Pscf