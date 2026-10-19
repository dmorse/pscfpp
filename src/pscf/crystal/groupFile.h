#ifndef PSCF_GROUP_FILE_H
#define PSCF_GROUP_FILE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
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
