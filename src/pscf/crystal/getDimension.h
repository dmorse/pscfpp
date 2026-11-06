#ifndef PSCF_GET_DIMENSION_H
#define PSCF_GET_DIMENSION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   /*
   * Get dimension of space from command line options.
   * 
   * \param argc number of command line arguments
   * \param argv vector of pointers to command line arguments
   * \return integer dimension of space
   */
   int getDimension(int argc, char **argv);

}
#endif
