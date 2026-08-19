/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Ensemble.h"

namespace Pscf { 

   using namespace Util;

   /* 
   * Extract a Ensemble from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, Ensemble& policy)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Closed" || buffer == "closed") {
         policy = Ensemble::Closed;
      } else 
      if (buffer == "Open" || buffer == "open") {
         policy = Ensemble::Open;
      } else {
         UTIL_THROW("Invalid Ensemble enum string in operator >>");
      } 
      return in;
   }
   
   /* 
   * Insert a Ensemble to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, Ensemble policy) 
   {
      if (policy == Ensemble::Closed) {
         out << "Closed";
      } else 
      if (policy == Ensemble::Open) {
         out << "Open";
      } else 
      if (policy == Ensemble::Unknown) {
         out << "Unknown";
      } else {
         std::cout << "Invalid Ensemble enum value on input" 
                   << std::endl;
         UTIL_THROW("Unrecognized value for Ensemble enum");
      } 
      return out;
   }

} // namespace Pscf
