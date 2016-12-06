/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include "Iterator.h"
#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
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
      domain_(),
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      p_(),
      m_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(0),
      hasDomain_(0),
      hasFields_(0)
   {  
      setClassName("System"); 

      #ifdef PSCF_GSL
      interactionPtr_ = new ChiInteraction(); 
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

      initHomogeneous();

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain());
      hasDomain_ = true;
      allocateFields();

      // Initialize iterator
      iterator().setSystem(*this);
      readParamComposite(in, iterator());
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
      UTIL_CHECK(hasDomain_);

      // Allocate memory in mixture
      mixture().setDomain(domain());

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      cFields_.allocate(nMonomer);
      int nx = domain().nx();
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
      //if (!isInitialized_) {
      //    UTIL_THROW("McSimulation is not initialized");
      //}

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
         if (command == "READ_WFIELDS") {
            inBuffer >> filename;
            Log::file() << Str(filename, 15) << std::endl;
            fileMaster().openInputFile(filename, inputFile);
            readWFields(inputFile);
            inputFile.close();
         } else
         if (command == "WRITE_WFIELDS") {
            inBuffer >> filename;
            Log::file() << Str(filename, 15) << std::endl;
            fileMaster().openOutputFile(filename, outputFile);
            writeFields(outputFile, wFields_);
            outputFile.close();
         } else
         if (command == "WRITE_CFIELDS") {
            inBuffer >> filename;
            Log::file() << Str(filename, 15) << std::endl;
            fileMaster().openOutputFile(filename, outputFile);
            writeFields(outputFile, cFields_);
            outputFile.close();
         } else
         if (command == "ITERATE") {
            Log::file() << std::endl;
            iterator().solve();
         } else 
         if (command == "COMPUTE_HOMOGENEOUS") {
            int mode;
            inBuffer >> mode;
            computeHomogeneous(mode);
            #if 0
            if (mode == 1) {
               double fHomo = homogeneous().fHelmholtz();
               double pHomo = homogeneous().pressure();
               double df = fHelmholtz() - fHomo;
               Log::file() << "f (homo) = " << fHomo << std::endl;
               Log::file() << "p (homo) = " << fHomo << std::endl;
               std::cout   << "Delta f    " << df    << std::endl;
            }
            #endif 
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

   /*
   * Compute Helmoltz free energy and pressure
   */
   void System::computeFreeEnergy()
   {
      fHelmholtz_ = 0.0;
 
      // Compute ideal gas contributions to fHelhmoltz_
      Polymer* polymerPtr;
      double phi, mu, length;
      int np = mixture().nPolymer();
      for (int i = 0; i < np; ++i) {
         polymerPtr = &mixture().polymer(i);
         phi = polymerPtr->phi();
         mu = polymerPtr->mu();
         // Recall: mu = ln(phi/q)
         length = polymerPtr->length();
         fHelmholtz_ += phi*( mu - 1.0 )/length;
      }

      // Apply Legendre transform subtraction
      int nm = mixture().nMonomer();
      for (int i = 0; i < nm; ++i) {
         fHelmholtz_ -= 
                  domain().innerProduct(wFields_[i], cFields_[i]);
      }

      // Add average interaction free energy density per monomer
      int nx = domain().nx();
      if (!f_.isAllocated()) f_.allocate(nx);
      if (!c_.isAllocated()) c_.allocate(nm);
      int j;
      for (int i = 0; i < nx; ++i) {
         for (j = 0; j < nm; ++j) {
            c_[j] = cFields_[j][i];
         }
         f_[i] = interaction().fHelmholtz(c_);
      }
      fHelmholtz_ = + domain().spatialAverage(f_);

      // Compute pressure
      pressure_ = -fHelmholtz_;
      for (int i = 0; i < np; ++i) {
         polymerPtr = & mixture().polymer(i);
         phi = polymerPtr->phi();
         mu = polymerPtr->mu();
         length = polymerPtr->length();
         pressure_ += phi*mu/length;
      }

   }

   void System::readWFields(std::istream &in)
   {
      UTIL_CHECK(hasDomain_);

      // Read grid dimensions
      std::string label;
      int nx, nm;
      in >> label;
      UTIL_CHECK(label == "nx");
      in >> nx;
      UTIL_CHECK(nx == domain().nx());
      in >> label;
      UTIL_CHECK (label == "nm");
      in >> nm;
      UTIL_CHECK(nm == mixture().nMonomer());

      // Read fields
      int i,j, idum;
      for (i = 0; i < nx; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nm; ++j) {
            in >> wFields_[j][i];
         }
      }
      double shift = wFields_[nm - 1][nx-1];
      for (i = 0; i < nx; ++i) {
         for (j = 0; j < nm; ++j) {
            wFields_[j][i] -= shift;
         }
      }

   }

   void System::writeFields(std::ostream &out, Array<Field> const &  fields)
   {
      int i, j;
      int nx = domain().nx();
      int nm = mixture().nMonomer();
      out << "nx     "  <<  nx              << std::endl;
      out << "nm     "  <<  nm              << std::endl;

      // Write fields
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][i]);
         }
         out << std::endl;
      }
   }

   void System::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      if (!c_.isAllocated()) c_.allocate(nm);
 
      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      for (i = 0; i < np; ++i) {
 
         // Count total number of monomers of each type
         for (j = 0; j < nm; ++j) {
            c_[j] = 0.0;
         }
         nb = mixture().polymer(i).nBlock(); 
         for (k = 0; k < nb; ++k) {
            Block& block = mixture().polymer(i).block(k);
            j = block.monomerId();
            c_[j] += block.length();
         }
 
         // Count the number of clumps of nonzero size
         nc = 0;
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 0.0) ++nc;
         }
         homogeneous_.molecule(i).setNClump(nc);
 
         // Set clump properties for this Homogeneous::Molecule
         k = 0; // Clump index
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 0.0) {
               homogeneous_.molecule(i).clump(k).setMonomerId(j);
               homogeneous_.molecule(i).clump(k).setSize(c_[j]);
               ++k;
            }
         }
        
      }

   }

   /*
   * Compute properties of a homogeneous reference system.
   *
   * mode == 0 : Composition equals spatial average composition
   * mode == 1:  Chemical potential equal to that of system,
   *             composition guess given at last grid point.
   * mode == 2:  Chemical potential equal to that of system,
   *             composition guess given at last grid point.
   */
   void System::computeHomogeneous(int mode)
   {
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      if (!p_.isAllocated()) p_.allocate(np+ns);
      if (mode == 0) {
         for (int i = 0; i < np; ++i) {
            p_[i] = mixture().polymer(i).phi(); 
         }
         homogeneous().setComposition(p_);
         double xi = 0.0;
         homogeneous().computeMu(interaction(), xi);
      } else 
      if (mode == 1 || mode == 2) {
         for (int i = 0; i < np; ++i) {
            m_[i] = mixture().polymer(i).mu(); 
         }
         int ix; // Grid index from which we obtain guess of 
         if (mode == 1) {
            ix = domain().nx() - 1;
         } else 
         if (mode == 2) {
            ix = 0;
         }
         for (int i = 0; i < np; ++i) {
            p_[i] = 0.0;
            int nb = mixture().polymer(i).nBlock();
            for (int j = 0; j < nb; ++j) {
               p_[i] += mixture().polymer(i).block(j).cField()[ix];
            }
         }
         double xi = 0.0;
         homogeneous().computePhi(interaction(), m_, p_, xi);
      } else {
         UTIL_THROW("Unknown mode in computeHomogeneous");
      }
   }

} // namespace Fd1d
} // namespace Pscf
