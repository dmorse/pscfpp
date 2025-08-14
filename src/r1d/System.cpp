/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <r1d/iterator/Iterator.h>
#include <r1d/iterator/IteratorFactory.h>
#include <r1d/sweep/Sweep.h>
#include <r1d/sweep/SweepFactory.h>
#include <r1d/iterator/NrIterator.h>
#include <r1d/misc/HomogeneousComparison.h>
#include <r1d/misc/FieldIo.h>

#include <pscf/inter/Interaction.h>
#include <pscf/inter/Interaction.h>
#include <pscf/floryHuggins/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/param/BracketPolicy.h>
#include <util/misc/ioUtil.h>

#include <string>
#include <unistd.h>

namespace Pscf {
namespace R1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   System::System()
    : mixture_(),
      domain_(),
      fileMaster_(),
      fieldIo_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      fIdeal_(0.0),
      fInter_(0.0),
      pressure_(0.0),
      hasMixture_(0),
      hasDomain_(0),
      hasFields_(0),
      hasSweep_(0)
   {  
      setClassName("System"); 

      fieldIo_.associate(domain_, fileMaster_);
      interactionPtr_ = new Interaction(); 
      iteratorFactoryPtr_ = new IteratorFactory(*this); 
      sweepFactoryPtr_ = new SweepFactory(*this);

      BracketPolicy::set(BracketPolicy::Optional);
   }

   /*
   * Destructor.
   */
   System::~System()
   {
      if (interactionPtr_) {
         delete interactionPtr_;
      }
      if (iteratorPtr_) {
         delete iteratorPtr_;
      }
      if (iteratorFactoryPtr_) {
         delete iteratorFactoryPtr_;
      }
      if (sweepPtr_) {
         delete sweepPtr_;
      }
      if (sweepFactoryPtr_) {
         delete sweepFactoryPtr_;
      }
   }

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
            oFlag = true;
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

      // Initialize FloryHuggins::Mixture object
      homogeneous_.initialize(mixture());

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain());
      hasDomain_ = true;
      allocateFields();

      // Instantiate and initialize an Iterator 
      std::string className;
      bool isEnd;
      iteratorPtr_ = iteratorFactoryPtr_->readObject(in, *this, 
                                                     className, isEnd);
      if (!iteratorPtr_) {
         std::string msg = "Unrecognized Iterator subclass name ";
         msg += className;
         UTIL_THROW(msg.c_str());
      }

      // Optionally instantiate and initialize a Sweep object
      sweepPtr_ = 
         sweepFactoryPtr_->readObjectOptional(in, *this, className, isEnd);
      if (sweepPtr_) {
         hasSweep_ = true;
         sweepPtr_->setSystem(*this);
      } else {
         hasSweep_ = false;
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
   * Write parameter file, omitting the sweep block.
   */
   void System::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture_.writeParam(out);
      interactionPtr_->writeParam(out);
      domain_.writeParam(out);
      iteratorPtr_->writeParam(out);
      out << "}" << std::endl;
   }

   /*
   * Allocate memory for fields.
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

      std::istream& inBuffer = in;

      bool readNext = true;
      while (readNext) {

         inBuffer >> command;

         if (inBuffer.eof()) {
            break;
         } else {
            Log::file() << command;
         }

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W") {
            readEcho(inBuffer, filename);
            fieldIo_.readFields(wFields(), filename);  
         } else
         if (command == "COMPUTE") {
            // Solve the modified diffusion equation, without iteration
            Log::file() << std::endl;
            compute();
         } else
         if (command == "ITERATE") {
            // Attempt iteratively solve a single SCFT problem
            bool isContinuation = false;
            int fail = iterate(isContinuation);
            if (fail) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            // Attempt to solve a sequence of SCFT problems along a line
            // through parameter space.
            sweep();
         } else
         if (command == "WRITE_PARAM") {
            readEcho(inBuffer, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file);
            writeParam(file);
            file.close();
         } else
         if (command == "WRITE_THERMO") {
            readEcho(inBuffer, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file, std::ios_base::app);
            writeThermo(file);
            file.close();
         } else
         if (command == "COMPARE_HOMOGENEOUS") {
            int mode;
            inBuffer >> mode;

            readEcho(inBuffer, filename);
            Log::file() << std::endl;
            Log::file() << "mode       = " << mode << std::endl;

            std::ofstream file;
            fileMaster().openOutputFile(filename, file, std::ios_base::app);

            HomogeneousComparison comparison(*this);
            comparison.compute(mode);
            comparison.output(mode, file);
            file.close();
         } else 
         if (command == "WRITE_W") {
            readEcho(inBuffer, filename);
            fieldIo_.writeFields(wFields(), filename);  
         } else
         if (command == "WRITE_C") {
            readEcho(inBuffer, filename);
            fieldIo_.writeFields(cFields(), filename);  
         } else
         if (command == "WRITE_BLOCK_C") {
            readEcho(inBuffer, filename);
            fieldIo_.writeBlockCFields(mixture_, filename);  
         } else
         if (command == "WRITE_Q_SLICE") {
            int polymerId, blockId, directionId, segmentId;
            readEcho(inBuffer, filename);
            inBuffer >> polymerId;
            inBuffer >> blockId;
            inBuffer >> directionId;
            inBuffer >> segmentId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n"
                        << Str("segment ID  ", 21) << segmentId << std::endl;
            writeQSlice(filename, polymerId, blockId, directionId, 
                                  segmentId);
         } else
         if (command == "WRITE_Q_TAIL") {
            readEcho(inBuffer, filename);
            int polymerId, blockId, directionId;
            inBuffer >> polymerId;
            inBuffer >> blockId;
            inBuffer >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQTail(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_VERTEX") {
            int polymerId, vertexId;
            inBuffer >> polymerId;
            Log::file() << std::endl;
            Log::file() << "polymerId = " 
                        << Int(polymerId, 5) << std::endl;
            inBuffer >> vertexId;
            Log::file() << "vertexId  = " 
                        << Int(vertexId, 5) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile   = " 
                        << Str(filename, 20) << std::endl;
            fieldIo_.writeVertexQ(mixture_, polymerId, vertexId, filename);  
         } else
         if (command == "WRITE_Q") {
            readEcho(inBuffer, filename);
            int polymerId, blockId, directionId;
            inBuffer >> polymerId;
            inBuffer >> blockId;
            inBuffer >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(inBuffer, filename);
            writeQAll(filename);
         } else
         if (command == "REMESH_W") {
            int nx;
            inBuffer >> nx;
            Log::file() << std::endl;
            Log::file() << "nx      = " << Int(nx, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo_.remesh(wFields(), nx, filename);
         } else
         if (command == "EXTEND_W") {
            int m;
            inBuffer >> m;
            Log::file() << std::endl;
            Log::file() << "m       = " << Int(m, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo_.extend(wFields(), m, filename);
         } else {
            Log::file() << "  Error: Unknown command  " 
                        << command << std::endl;
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
   * Read a string (e.g., a filename) and echo to the log file.
   */
   void System::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      Log::file() << "  " << Str(string, 20) << std::endl;
   }

   // Primary SCFT Computations

   /*
   * Solve MDE for current w-fields, without iteration.
   */
   void System::compute()
   {
      //UTIL_CHECK(hasWFields_);

      // Solve the modified diffusion equation (without iteration)
      mixture_.compute(wFields(), cFields_);

      //hasCFields_ = true;
   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   int System::iterate(bool isContinuation)
   {
      // UTIL_CHECK(hasWFields_);
      // hasCFields_ = false;

      Log::file() << std::endl;

      // Call iterator (return 0 for convergence, 1 for failure)
      int error = iterator().solve(isContinuation);
      
      //hasCFields_ = true;

      // If converged, compute related properties
      if (!error) {   
         computeFreeEnergy();
         writeThermo(Log::file());
      }
      return error;
   }

   /*
   * Perform sweep along a line in parameter space.
   */
   void System::sweep()
   {
      //UTIL_CHECK(hasWFields_);
      UTIL_CHECK(sweepPtr_);
      Log::file() << std::endl;
      Log::file() << std::endl;

      sweepPtr_->sweep();
   }
  
   // Thermodynamic properties
 
   /*
   * Compute Helmholtz free energy and pressure
   */
   void System::computeFreeEnergy()
   {
      fHelmholtz_ = 0.0;
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      int nm = mixture().nMonomer();
      double phi, mu;
 
      // Polymeric ideal gas contributions to fHelhmoltz_
      if (np > 0) {
         Polymer* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture().polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            // Recall: mu = ln(phi/q)
            length = polymerPtr->length();
            if (phi > 1E-08) {
               fHelmholtz_ += phi*( mu - 1.0 )/length;
            }
         }
      }

      // Solvent ideal gas contributions to fHelhmoltz_
      if (ns > 0) {
         Solvent* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture().solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            // Recall: mu = ln(phi/q)
            size = solventPtr->size();
            if (phi > 1E-08) {
               fHelmholtz_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      // Apply Legendre transform subtraction
      for (int i = 0; i < nm; ++i) {
         fHelmholtz_ -= 
                  domain().innerProduct(wFields_[i], cFields_[i]);
      }
      fIdeal_ = fHelmholtz_;

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
      fInter_ = fHelmholtz_ - fIdeal_;

      // Compute pressure
      pressure_ = -fHelmholtz_;
      if (np > 0) {
         double length;
         Polymer* polymerPtr;
         for (int i = 0; i < np; ++i) {
            polymerPtr = & mixture().polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            if (phi > 1E-08) {
               pressure_ += phi*mu/length;
            }
         }
      }
      if (ns > 0) {
         double size;
         Solvent* solventPtr;
         for (int i = 0; i < ns; ++i) {
            solventPtr = & mixture().solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1E-08) {
               pressure_ += phi*mu/size;
            }
         }
      }

   }

   void System::writeThermo(std::ostream& out) const
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
      out << std::endl;

      // Polymers
      int np = mixture().nPolymer();
      if (np > 0) {
         out << "polymers:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < np; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
                << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

      // Solvents
      int ns = mixture().nSolvent();
      if (ns > 0) {
         out << "solvents:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < ns; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture().solvent(i).phi(),18, 11)
                << "  " << Dbl(mixture().solvent(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

   }

   // Output operations (correspond to command file commands)

   /*
   * Write w-fields in symmetry-adapted basis format. 
   */
   void System::writeW(std::string const & filename) const
   {
      //UTIL_CHECK(hasWFields_);
      fieldIo_.writeFields(wFields(), filename);
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   void System::writeC(std::string const & filename) const
   {
      //UTIL_CHECK(hasCFields_);
      fieldIo_.writeFields(cFields(), filename);
   }

   /*
   * Write all concentration fields in real space (r-grid) format, for 
   * each block (or solvent) individually rather than for each species.
   */
   void System::writeBlockC(std::string const & filename) const
   {
      //UTIL_CHECK(hasCFields_);
      fieldIo_.writeBlockCFields(mixture_, filename);
   }
   
   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   void System::writeQSlice(const std::string & filename, 
                            int polymerId, int blockId, 
                            int directionId, int segmentId) 
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator const & 
           propagator = polymer.propagator(blockId, directionId);
      Field const& field = propagator.q(segmentId);
      fieldIo_.writeField(field, filename);
   }

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   void System::writeQTail(const std::string & filename, 
                           int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Field const& 
            field = polymer.propagator(blockId, directionId).tail();
      fieldIo_.writeField(field, filename);
   }

   /*
   * Write the propagator for a block and direction.
   */
   void System::writeQ(const std::string & filename, 
                       int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator const& propagator 
                              = polymer.propagator(blockId, directionId);
      int nslice = propagator.ns();
      
      // Open file
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);

      // Write header
      file << "ngrid" << std::endl
           << "          " << domain_.nx() << std::endl
           << "nslice"    << std::endl
           << "          " << nslice << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < nslice; ++i) {
          file << "slice " << i << std::endl;
          Field const& field = propagator.q(i);
          fieldIo_.writeField(field, file, hasHeader);
      }
      file.close();
   }

   /*
   * Write propagators for all blocks of all polymers to files.
   */
   void System::writeQAll(std::string const & basename) const
   {
      std::string filename;
      int np, nb, ip, ib, id;
      np = mixture_.nPolymer();
      for (ip = 0; ip < np; ++ip) {
         nb = mixture_.polymer(ip).nBlock();
         for (ib = 0; ib < nb; ++ib) {
            for (id = 0; id < 2; ++id) {
               filename = basename;
               filename += "_";
               filename += toString(ip);
               filename += "_";
               filename += toString(ib);
               filename += "_";
               filename += toString(id);
               //filename += ".q";
               writeQ(filename, ip, ib, id);
            }
         }
      }
   }

} // namespace R1d
} // namespace Pscf
