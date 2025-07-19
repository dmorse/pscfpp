/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/system/System.h>
#include <prdc/crystal/getDimension.h>
#include <pscf/chem/PolymerModel.h>

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

      // Prohibit later changes to global polymer model
      PolymerModel::lock();

      // Read command script and execute commands
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
