/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "groupFile.h"

#define AS_STRING(s) # s

namespace Pscf { 
namespace Pssp 
{ 

   /**
   * Generates the file name from a group name.
   *
   * \param groupName standard name of space group
   */
   std::string makeGroupFileName(int D, std::string groupName)
   {
      std::string filename = AS_STRING(DAT_DIR) ;
      filename += "/groups/";
      if (D==1) {
         filename += "1/";
      } else
      if (D==2) {
         filename += "2/";
      } else
      if (D==2) {
         filename += "3/";
      }
      filename += groupName;
      return filename;
   }

} // namespace Pscf:Pssp
} // namespace Pscf
