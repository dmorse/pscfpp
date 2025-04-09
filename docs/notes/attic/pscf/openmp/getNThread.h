#ifndef PRDC_GET_N_THREAD_H
#define PRDC_GET_N_THREAD_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf{

   /*
   * Get number of threads from command line option -t and set default.
   * 
   * This function searches the argument list argv for option -t and,
   * if found, sets the default number of OpenMP threads to the argument
   * of this option.
   * 
   * \param argc number of command line arguments
   * \param argv vector of pointers to command line arguments
   * \return integer dimension of space
   */
   int getNThread(int argc, char **argv);

}
#endif
