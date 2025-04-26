/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParameterModifier.h"

namespace Pscf {

   // Constructor
   ParameterModifier::ParameterModifier()
   {}

   // Destructor
   ParameterModifier::~ParameterModifier()
   {}

   /*
   * Return specialized sweep parameter types to add to the Sweep object.
   */
   GArray<ParameterType> ParameterModifier::getParameterTypes()
   {  
      // Default behavior: return empty array
      GArray<ParameterType> arr;
      return arr; 
   } 

   /*
   * Set the value of a specialized sweep parameter, and throw an error 
   * if the parameter is not found.
   */
   void ParameterModifier::setParameter(std::string name, 
                                        DArray<int> ids, double value)
   {
      bool success(true);
      setParameter(name, ids, value, success);
      if (!success) {
         UTIL_THROW(("Parameter name " + name + " not recognized.").c_str());
      }
   }

   /*
   * Get the value of a specialized sweep parameter, and throw an error 
   * if the parameter is not found.
   */
   double ParameterModifier::getParameter(std::string name, DArray<int> ids)
   const
   {
      bool success(true);
      double val = getParameter(name, ids, success);
      if (!success) {
         UTIL_THROW(("Parameter name " + name + " not recognized.").c_str());
      }
      return val;
   }

}