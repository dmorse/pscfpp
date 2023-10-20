#ifndef PRDC_GET_N_THREAD_H
#define PRDC_GET_N_THREAD_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf{
namespace Prdc{

   /*
   * Get the number of threads from command line option -t.
   * 
   * \param argc number of command line arguments
   * \param argv vector of pointers to command line arguments
   * \return integer dimension of space
   */
   int getNThread(int argc, char **argv);

}
}
#endif
