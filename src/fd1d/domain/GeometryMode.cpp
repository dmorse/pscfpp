/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>     // uses UTIL_THROW
#include "GeometryMode.h"    // class header

namespace Pscf{
namespace Fd1d
{

   using namespace Util;

   /* 
   * Extract a GeometryMode from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, GeometryMode& lattice)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Planar" || buffer == "planar") {
         lattice = Planar;
      } else 
      if (buffer == "Cylindrical" || buffer == "cylindrical") {
         lattice = Cylindrical;
      } else 
      if (buffer == "Spherical" || buffer == "spherical") {
         lattice = Spherical;
      } else {
         UTIL_THROW("Invalid GeometryMode value input");
      }
      return in;
   }
   
   /* 
   * Insert a GeometryMode to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, GeometryMode lattice) 
   {
      if (lattice == Planar) {
         out << "planar";
      } else 
      if (lattice == Cylindrical) {
         out << "cylindrical";
      } else
      if (lattice == Spherical) {
         out << "spherical";
      } else {
         UTIL_THROW("This should never happen");
      } 
      return out; 
   }

}
}
