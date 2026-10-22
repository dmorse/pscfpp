/// \file 

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <r1/system/System.h>

/**
* Main function for pscf_1d program.
*
* \ingroup Pscf_R1_Module
*
* \param argc number of command line parameter strings
* \param argv array of command line parameter strings
*/
int main(int argc, char **argv)
{
   Pscf::R1::System system;

   // Process command line options
   system.setOptions(argc, argv);

   // Read parameters from default parameter file
   system.readParam();

   // Read command script to run system
   system.readCommands();

   return 0;
}
