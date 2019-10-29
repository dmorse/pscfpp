/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fd1d/System.h>

/**
* \page pscf_fd1d_page pscf_fd1d
*
* Polymer Self-Consistent Field Theory - One-dimensional Finite Difference
*
* Usage:
*
*    pscf_fd1d [-e] [-r file] [-p file] [-c file] [-i prefix] [-o prefix]
*
* Options:
*
*   -e  
*
*    Enable echoing of parameter file to log file as it is read. This
*    option is often useful for debugging the parameter file.
*
*  -p file
*
*   Set the parameter file name, given by the argument "file". 
*   The -p and -r options are incompatible. 
*
*  -c file
*
*   Set the command file name, given by the argument "file". The
*   command file may also be specified in the parameter file.
*
*  -i prefix
*
*   Set the input file path prefix, given by the argument "prefix".
*
*  -o prefix
*
*   Set the output file path prefix, given by the argument "prefix".
*
*   Set replicated mode for parallel simulations.
*
*/

int main(int argc, char **argv)
{
   Pscf::Fd1d::System system;

   // Process command line options
   system.setOptions(argc, argv);

   // Read parameters from default parameter file
   system.readParam();

   // Read command script to run system
   system.readCommands();

   return 0;
}
