#ifndef PRDC_GROUP_FILE_H
#define PRDC_GROUP_FILE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

namespace Pscf { 
namespace Prdc { 

   /**
   * Generates the file name from a group name.
   *
   * \param D dimensionality of space (D=1,2 or 3)
   * \param groupName standard name of space group
   * \ingroup Prdc_Crystal_Module
   */
   std::string makeGroupFileName(int D, std::string groupName);


} // namespace Prdc
} // namespace Pscf
#endif
