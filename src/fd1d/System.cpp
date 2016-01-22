/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include "Iterator.h"
#include "NrIterator.h"

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   System::System()
    : mixture_(),
      iteratorPtr_(0)
   {  
      setClassName("System"); 
      iteratorPtr_ = new NrIterator(); 
   }

   /*
   * Process command line options.
   */
   void System::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      #ifdef MCMD_PERTURB
      bool  fflag = false;  // free energy perturbation
      #endif
      char* pArg = 0;
      char* cArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f")) != -1) {
         switch (c) {
         case 'e':
            eflag = true;
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }
   
      // Set flag to echo parameters as they are read.
      if (eflag) {
         Util::ParamComponent::setEcho(true);
      }

      // If option -p, set parameter file name
      if (pFlag) {
         fileMaster().setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      }

      #if 0
      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

      // If option -r, restart
      if (rFlag) {
         //Log::file() << "Reading restart file " 
         //            << std::string(rarg) << std::endl;
         isRestarting_ = true; 
         load(std::string(rarg));
      }
      #endif

   }

   /*
   * Read parameters and initialize.
   */
   void System::readParameters(std::istream& in)
   {
   }  

} // namespace Fd1d
} // namespace Pscf
