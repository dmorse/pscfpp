/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <fd1d/iterator/Iterator.h>
#include <fd1d/sweep/Sweep.h>
#include <fd1d/sweep/SweepFactory.h>
#ifdef PSCF_GSL
#include <fd1d/iterator/NrIterator.h>
#endif
#include <fd1d/misc/HomogeneousComparison.h>
#include <fd1d/misc/FieldIo.h>

#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

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
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(0),
      hasDomain_(0),
      hasFields_(0),
      hasSweep_(0)
   {  
      setClassName("System"); 

      #ifdef PSCF_GSL
      interactionPtr_ = new ChiInteraction(); 
      iteratorPtr_ = new NrIterator(*this); 
      #endif
      sweepFactoryPtr_ = new SweepFactory(*this);
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

      // Initialize homogeneous object
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain());
      hasDomain_ = true;
      allocateFields();

      // Initialize iterator
      readParamComposite(in, iterator());

      // Optionally instantiate a Sweep object
      readOptional<bool>(in, "hasSweep", hasSweep_);
      if (hasSweep_) {
         std::string className;
         bool isEnd;
         sweepPtr_ = 
            sweepFactoryPtr_->readObject(in, *this, className, isEnd);
         if (!sweepPtr_) {
            UTIL_THROW("Unrecognized Sweep subclass name");
         }
         sweepPtr_->setSystem(*this);
      }
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
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasDomain_);

      std::string command;
      std::string filename;
      FieldIo fieldIo(*this);

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
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo.readFields(filename, wFields());  
         } else
         if (command == "ITERATE") {
            Log::file() << std::endl;
            Log::file() << std::endl;
            iterator().solve();
            outputThermo(Log::file());
         } else 
         if (command == "COMPARE_HOMOGENEOUS") {
            int mode;
            inBuffer >> mode;
            Log::file() << std::endl;
            Log::file() << "mode       = " << mode << std::endl;

            HomogeneousComparison comparison(*this);
            comparison.compute(mode);
            comparison.output(mode, Log::file());
         } else 
         if (command == "SWEEP") {
            UTIL_CHECK(hasSweep_);
            UTIL_CHECK(sweepPtr_);
            sweepPtr_->solve();
         } else 
         if (command == "WRITE_WFIELDS") {
            inBuffer >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo.writeFields(filename, wFields());  
         } else
         if (command == "WRITE_CFIELDS") {
            inBuffer >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo.writeFields(filename, cFields());  
         } else
         if (command == "REMESH_WFIELDS") {
            int nx;
            inBuffer >> nx;
            Log::file() << std::endl;
            Log::file() << "nx      = " << Int(nx, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo.remesh(filename, wFields(), nx);
         } else
         if (command == "EXTEND_WFIELDS") {
            int m;
            inBuffer >> m;
            Log::file() << std::endl;
            Log::file() << "m       = " << Int(m, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo.extend(filename, wFields(), m);
         } else
         {
            Log::file() << "  Error: Unknown command  " << command << std::endl;
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
         // Store c_[j] = local concentration of species j
         for (j = 0; j < nm; ++j) {
            c_[j] = cFields_[j][i];
         }
         // Compute f_[i] = excess free eenrgy at grid point i
         f_[i] = interaction().fHelmholtz(c_);
      }
      fHelmholtz_ += domain().spatialAverage(f_);

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

   #if 0
   void System::readWFields(std::istream &in)
   {
      UTIL_CHECK(hasDomain_);

      // Read grid dimensions
      std::string label;
      int nx, nm;
      in >> label;
      UTIL_CHECK(label == "nx");
      in >> nx;
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(nx == domain().nx());
      in >> label;
      UTIL_CHECK (label == "nm");
      in >> nm;
      UTIL_CHECK(nm > 0);
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

   }

   void System::writeFields(std::ostream &out, 
                            Array<Field> const &  fields)
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
            out << "  " << Dbl(fields[j][i], 18, 11);
         }
         out << std::endl;
      }
   }
   #endif

   void System::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      // Allocate c_ work array, if necessary
      if (c_.isAllocated()) {
         UTIL_CHECK(c_.capacity() == nm);
      } else {
         c_.allocate(nm);
      }

      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      for (i = 0; i < np; ++i) {

         // Initial array of clump sizes 
         for (j = 0; j < nm; ++j) {
            c_[j] = 0.0;
         }

         // Compute clump sizes for all monomer types.
         nb = mixture().polymer(i).nBlock(); 
         for (k = 0; k < nb; ++k) {
            Block& block = mixture().polymer(i).block(k);
            j = block.monomerId();
            c_[j] += block.length();
         }
 
         // Count the number of clumps of nonzero size
         nc = 0;
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               ++nc;
            }
         }
         homogeneous_.molecule(i).setNClump(nc);
 
         // Set clump properties for this Homogeneous::Molecule
         k = 0; // Clump index
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               homogeneous_.molecule(i).clump(k).setMonomerId(j);
               homogeneous_.molecule(i).clump(k).setSize(c_[j]);
               ++k;
            }
         }
         homogeneous_.molecule(i).computeSize();

      }

   }

   void System::outputThermo(std::ostream& out)
   {
      out << std::endl;
      out << "fHelmholtz = " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure   = " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;

      out << "Polymers:" << std::endl;
      out << "    i"
          << "        phi[i]      "
          << "        mu[i]       " 
          << std::endl;
      for (int i = 0; i < mixture().nPolymer(); ++i) {
         out << Int(i, 5) 
             << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
             << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)  
             << std::endl;
      }
      out << std::endl;
   }

} // namespace Fd1d
} // namespace Pscf
