/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/getDimension.h>
#include <rpc/System.h>

#include <iostream>

namespace Pscf {
namespace Rpc {

   /**
   * Function template for main pscf_pc program.
   *
   * \param argc  number of command line parameters
   * \param argv  array of command line parameter strings
   * \ingroup Pscf_Rpc_Module
   */
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

/**
* Main pscf_pc program.
*
* \param argc  number of command line arguments
* \param argv  array of command line arguments
* \ingroup Pscf_Rpc_Module
*/
int main(int argc, char **argv)
{

   // Extract the dimension of space from argument of -d option
   int D = Pscf::Prdc::getDimension(argc, argv);
   std::cout << "dimension   " << D << std::endl;

   if (1 == D) {
      Pscf::Rpc::run<1>(argc, argv);
   } else
   if (2 == D) {
      Pscf::Rpc::run<2>(argc, argv);
   } else
   if (3 == D) {
      Pscf::Rpc::run<3>(argc, argv);
   } else {
      std::cout << " Invalid dimension = " << D << std::endl;
   }

   fftw_cleanup();
}
