#ifndef PSSP_GROUP_FILE_H
#define PSSP_GROUP_FILE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

namespace Pscf { 

   /**
   * Generates the file name from a group name.
   *
   * \param D dimensionality of space (D=1,2 or 3)
   * \param groupName standard name of space group
   * \ingroup Pscf_Crystal_Module
   */
   std::string makeGroupFileName(int D, std::string groupName);


} // namespace Pscf
#endif
