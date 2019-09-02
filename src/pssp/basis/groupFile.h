#ifndef PSSP_GROUP_FILE_H
#define PSSP_GROUP_FILE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

namespace Pscf { 
namespace Pssp 
{ 

   /**
   * Generates the file name from a group name.
   *
   * \param groupName standard name of space group
   */
   std::string makeGroupFileName(int D, std::string groupName);


} // namespace Pscf::Pssp
} // namespace Pscf
#endif
