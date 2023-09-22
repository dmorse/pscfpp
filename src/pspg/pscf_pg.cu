/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/mesh/getDimension.h>
#include <pspg/System.h>
#include <iostream>

namespace Pscf {
namespace Pspg {

   template <int D>
   void run(int argc, char **argv) {
      System<D> system;

      // Process command line options
      system.setOptions(argc, argv);

      // Read parameters from default parameter file
      system.readParam();

      // Read command script to run system
      system.readCommands();
   }

}
}

/// Main program
int main(int argc, char **argv)
{

   // Extract the dimension of space from argument of -d option
   int D;
   char** argvcopy = argv;
   D = Pscf::getDimension(argc, argvcopy);
   std::cout << "dimension    " << D << std::endl;

   if (1 == D) {
      Pscf::Pspg::run<1>(argc, argv);
   } else
   if (2 == D) {
      Pscf::Pspg::run<2>(argc, argv);
   } else
   if (3 == D) {
      Pscf::Pspg::run<3>(argc, argv);
   } else {
      std::cout << " Invalid dimension = " << D << std::endl;
   }

}
