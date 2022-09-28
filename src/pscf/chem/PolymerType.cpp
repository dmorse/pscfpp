/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerType.h"
#include <util/global.h>
#include <string>

namespace Pscf
{ 

   using namespace Util;

   /*
   * Input stream extractor for a PolymerType enumeration.
   */ 
   std::istream& operator >> (std::istream& in, PolymerType::Enum& type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Branched" || buffer == "branched") {
         type = PolymerType::Branched;
      } else
      if (buffer == "Linear" || buffer == "linear") {
         type = PolymerType::Linear;
      } else {
         std::string msg = "Unknown input PolymerType value string: ";
         msg += buffer;
         UTIL_THROW(msg.c_str());
      } 
      return in;
   }

   /*
   * Input stream extractor for a PolymerType enumeration.
   */ 
   std::ostream& operator << (std::ostream& out, PolymerType::Enum& type)
   {
      if (type == PolymerType::Branched) {
         out << "branched";
      } else
      if (type == PolymerType::Linear) {
         out << "linear";
      } else {
         // This should never happen
         UTIL_THROW("Error writing a PolymerType value");
      }
      return out;
   }

}
