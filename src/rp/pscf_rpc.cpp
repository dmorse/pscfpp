/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/system/System.h>
#include <pscf/cpu/CppTp.h>
#include <prdc/crystal/getDimension.h>
#include <pscf/chem/PolymerModel.h>
#include <util/format/Int.h>
#include <util/misc/Log.h>
#include <util/misc/Memory.h>

#include <iostream>

namespace Pscf {
namespace Rpc {

   /**
   * Function template for main pscf_rpc program.
   *
   * \ingroup Pscf_Rpc_Module
   *
   * \param argc  number of command line parameters
   * \param argv  array of command line parameter strings
   */
   template <int D>
   void run(int argc, char **argv) {
      Rp::System<D, CppTp<D> > system;

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
* Main pscf_rpc program.
*
* \param argc  number of command line arguments
* \param argv  array of command line arguments
* \ingroup Pscf_Rpc_Module
*/
int main(int argc, char **argv)
{

   using namespace Util;
   using namespace Pscf;

   // Extract the dimension of space from argument of -d option
   int D = Pscf::Prdc::getDimension(argc, argv);
   std::cout << "dimension   " << D << std::endl;

   // Run relevant template specialization
   if (1 == D) {
      Rpc::run<1>(argc, argv);
   } else
   if (2 == D) {
      Rpc::run<2>(argc, argv);
   } else
   if (3 == D) {
      Rpc::run<3>(argc, argv);
   } else {
      std::cout << " Invalid dimension = " << D << std::endl;
   }

   // Report memory usage
   Log::file() << "\nCPU heap array memory usage:";
   Log::file() << "\n  Maximum array memory = " 
               << Int(Memory::max(), 12) << " bytes";
   Log::file() << "\n  Final array memory   = " 
               << Int(Memory::total(), 12) << " bytes";
   Log::file() << std::endl;

   fftw_cleanup();
}
