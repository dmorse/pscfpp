/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include "Iterator.h"
#include <util/format/Str.h>
#ifdef PSCF_GSL
#include "NrIterator.h"
#endif

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
      grid_(),
      fileMaster_(),
      iteratorPtr_(0),
      hasMixture_(0),
      hasGrid_(0),
      hasFields_(0)
   {  
      setClassName("System"); 

      #ifdef PSCF_GSL
      iteratorPtr_ = new NrIterator(); 
      #endif
   }

   /*
   * Destructor.
   */
   System::~System()
   {}

   /*
   * Process command line options.
   */
   void System::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
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
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            iFlag = true;
            oArg  = optarg;
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

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

   }

   /*
   * Read parameters and initialize.
   */
   void System::readParameters(std::istream& in)
   {
      readParamComposite(in, mixture());
      hasMixture_ = true;

      readParamComposite(in, grid());
      hasGrid_ = true;
      allocateFields();

      // Initialize iterator
      iterator().setSystem(*this);
   }

   /*
   * Read default parameter file.
   */
   void System::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());  
      readParameters(in);  
      readEnd(in);  
   }

   /*
   * Read default parameter file.
   */
   void System::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read parameters and initialize.
   */
   void System::allocateFields()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasGrid_);
      UTIL_CHECK(!hasFields_);

      mixture().setGrid(grid());
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      cFields_.allocate(nMonomer);
      int nx = grid().nx();
      for (int i = 0; i < nMonomer; ++i) {
         wField(i).allocate(nx);
         cField(i).allocate(nx);
      }
      hasFields_ = true;
   }

   /*
   * Read and execute commands from a specified command file.
   */
   void System::readCommands(std::istream &in)
   {
      // if (!isInitialized_) {
      //    UTIL_THROW("McSimulation is not initialized");
      // }

      std::string command;
      std::string filename;
      std::ifstream inputFile;
      std::ofstream outputFile;

      std::istream& inBuffer = in;

      bool readNext = true;
      while (readNext) {

         inBuffer >> command;
         Log::file() << command;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_OMEGA") {
            inBuffer >> filename;
            Log::file() << Str(filename, 15) << std::endl;
            fileMaster().openInputFile(filename, inputFile);
            readOmega(inputFile);
            inputFile.close();
         } else
         if (command == "ITERATE") {
            iterator().solve();
         } else 
         {
            Log::file() << "  Error: Unknown command  " << std::endl;
            readNext = false;
         }

      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   void System::readCommands()
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }

   void System::readOmega(std::istream &in)
   {

      // Read grid parameters:
      std::string label;
      double xMin, xMax;
      int nx;
      in >> label;
      UTIL_CHECK (label != "xMin");
      in >> xMin;
      in >> label;
      UTIL_CHECK (label != "xMax");
      in >> xMax;
      in >> label;
      UTIL_CHECK (label != "nx");
      in >> nx;
      if (!hasGrid_) {
         grid().setParameters(xMin, xMax, nx);
         allocateFields();
      } else {
         UTIL_CHECK(nx == grid().nx());
      }

      // Read fields
      int i,j, idum;
      int nMonomer = mixture().nMonomer();
      for (i = 0; i < nx; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nMonomer; ++j) {
            in >> wFields_[i][j];
         }
      }

   }

   void System::writeOmega(std::ostream &out)
   {
   }

} // namespace Fd1d
} // namespace Pscf
