#ifndef PRDC_GET_DIMENSION_H
#define PRDC_GET_DIMENSION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf{
namespace Prdc{

   /*
   * Get dimension of space from command line options.
   * 
   * \param argc number of command line arguments
   * \param argv vector of pointers to command line arguments
   * \return integer dimension of space
   */
   int getDimension(int argc, char **argv);

}
}
#endif
