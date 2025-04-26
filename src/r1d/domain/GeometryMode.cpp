/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>     // uses UTIL_THROW
#include "GeometryMode.h"    // class header

namespace Pscf{
namespace R1d
{

   using namespace Util;

   /* 
   * Extract a GeometryMode from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, GeometryMode& mode)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Planar" || buffer == "planar") {
         mode = Planar;
      } else 
      if (buffer == "Cylindrical" || buffer == "cylindrical") {
         mode = Cylindrical;
      } else 
      if (buffer == "Spherical" || buffer == "spherical") {
         mode = Spherical;
      } else {
         UTIL_THROW("Invalid GeometryMode value input");
      }
      return in;
   }
   
   /* 
   * Insert a GeometryMode to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, GeometryMode mode) 
   {
      if (mode == Planar) {
         out << "planar";
      } else 
      if (mode == Cylindrical) {
         out << "cylindrical";
      } else
      if (mode == Spherical) {
         out << "spherical";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out; 
   }

}
}
