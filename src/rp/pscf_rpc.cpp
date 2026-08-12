/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/run.h>
#include <prdc/crystal/getDimension.h>
#include <pscf/backends/CPT.h>
#include <iostream>

/**
* Main pscf_rpg program.
*
* \param argc  number of command line arguments
* \param argv  array of command line arguments
* \ingroup Pscf_Rp_Module
*/
int main(int argc, char **argv)
{

   // Extract the dimension of space from argument of -d option
   int D = Pscf::Prdc::getDimension(argc, argv);
   std::cout << "dimension    " << D << std::endl;

   if (1 == D) {
      Pscf::Rp::run<1, Pscf::CPT>(argc, argv);
   } else
   if (2 == D) {
      Pscf::Rp::run<2, Pscf::CPT>(argc, argv);
   } else
   if (3 == D) {
      Pscf::Rp::run<3, Pscf::CPT>(argc, argv);
   } else {
      std::cout << " Invalid dimension = " << D << std::endl;
   }

}
