/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.h>

namespace Pscf {
namespace Rp {

   /**
   * Function template for main pscf_rp* programs.
   *
   * \ingroup Pscf_Rp_Module
   *
   * \param argc  number of command line parameters
   * \param argv  array of command line parameter strings
   */
   template <int D, class T>
   void run(int argc, char **argv) {

      // Construct System object for specific dimension D.
      System<D,T> system;

      // Process command line options
      system.setOptions(argc, argv);

      // Read parameters from default parameter file
      system.readParam();

      // Read command script to run system
      system.readCommands();

   }

}
}

