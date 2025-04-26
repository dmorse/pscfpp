/// \file 

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <r1d/System.h>

/**
* Main function for pscf_1d program.
*
* \param argc number of command line parameter strings
* \param argv array of command line parameter strings
* \ingroup Pscf_R1d_Module
*/
int main(int argc, char **argv)
{
   Pscf::R1d::System system;

   // Process command line options
   system.setOptions(argc, argv);

   // Read parameters from default parameter file
   system.readParam();

   // Read command script to run system
   system.readCommands();

   return 0;
}
