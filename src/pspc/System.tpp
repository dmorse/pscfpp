#ifndef PSPC_SYSTEM_TPP
#define PSPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <pspc/sweep/Sweep.h>
#include <pspc/compressor/Compressor.h>

#include <pspc/iterator/IteratorFactory.h>
#include <pspc/sweep/SweepFactory.h>
#include <pspc/compressor/CompressorFactory.h>

#include <pspc/solvers/Mixture.h>
#include <pspc/solvers/Polymer.h>
#include <pspc/solvers/Solvent.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/inter/Interaction.h>
#include <pscf/inter/Interaction.h>
#include <pscf/homogeneous/Clump.h>

#include <pspc/field/BFieldComparison.h>
#include <pspc/field/RFieldComparison.h>

#include <util/param/BracketPolicy.h>
#include <util/misc/ioUtil.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   System<D>::System()
    : mixture_(),
      domain_(),
      mcMoveManager_(*this),
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      compressorPtr_(0),
      compressorFactoryPtr_(0),
      w_(),
      c_(),
      h_(),
      mask_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocatedRGrid_(false),
      isAllocatedBasis_(false)
   {  
      setClassName("System"); 
      domain_.setFileMaster(fileMaster_);
      w_.setFieldIo(domain_.fieldIo());
      h_.setFieldIo(domain_.fieldIo());
      mask_.setFieldIo(domain_.fieldIo());
      interactionPtr_ = new Interaction(); 
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this); 
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      compressorFactoryPtr_ = new CompressorFactory<D>(*this);
      BracketPolicy::set(BracketPolicy::Optional);
   }

   /*
   * Destructor.
   */
   template <int D>
   System<D>::~System()
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
      if (compressorPtr_) {
         delete compressorPtr_;
      }
      if (compressorFactoryPtr_) {
         delete compressorFactoryPtr_;
      }
   }

   /*
   * Process command line options.
   */
   template <int D>
   void System<D>::setOptions(int argc, char **argv)
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
         fileMaster_.setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster_.setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster_.setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster_.setOutputPrefix(std::string(oArg));
      }

   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {
      // Read the Mixture{ ... } block
      readParamComposite(in, mixture_);
      hasMixture_ = true;

      int nm = mixture_.nMonomer(); 
      int np = mixture_.nPolymer(); 
      int ns = mixture_.nSolvent(); 
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(np >= 0);
      UTIL_CHECK(ns >= 0);
      UTIL_CHECK(np + ns > 0);

      // Initialize homogeneous object 
      // NOTE: THIS OBJECT IS NOT USED AT ALL.
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      // Read the Interaction{ ... } block
      interaction().setNMonomer(nm);
      readParamComposite(in, interaction());

      // Read the Domain{ ... } block
      readParamComposite(in, domain_);

      mixture_.setMesh(domain_.mesh());
      mixture_.setupUnitCell(unitCell());

      // Allocate field array members of System
      allocateFieldsGrid();
      if (domain_.basis().isInitialized()) {
         allocateFieldsBasis();
      }

      // Optionally instantiate an Iterator object
      std::string className;
      bool isEnd;
      iteratorPtr_ = 
         iteratorFactoryPtr_->readObjectOptional(in, *this, className, 
                                                 isEnd);
      if (!iteratorPtr_) {
         Log::file() << "Notification: No iterator was constructed\n";
      }

      // Optionally instantiate a Sweep object
      if (iteratorPtr_) {
         sweepPtr_ = 
            sweepFactoryPtr_->readObjectOptional(in, *this, className, 
                                                 isEnd);
      }

      // Optionally instantiate a Compressor object
      compressorPtr_ = 
         compressorFactoryPtr_->readObjectOptional(in, *this, className, isEnd);

   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());
      readParameters(in);  
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam()
   {  readParam(fileMaster_.paramFile()); }

   /*
   * Write parameter file, omitting any sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture().writeParam(out);
      interaction().writeParam(out);
      domain_.writeParam(out);
      if (iteratorPtr_) {
         iterator().writeParam(out);
      }
      out << "}" << std::endl;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocateFieldsGrid()
   {
      // Preconditions 
      UTIL_CHECK(hasMixture_);
      int nMonomer = mixture_.nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(!isAllocatedRGrid_);

      // Alias for mesh dimensions
      IntVec<D> const & dimensions = domain_.mesh().dimensions();

      // Allocate w (chemical potential) fields
      w_.setNMonomer(nMonomer);
      w_.allocateRGrid(dimensions);

      // Allocate c (monomer concentration) fields
      c_.setNMonomer(nMonomer);
      c_.allocateRGrid(dimensions);

      h_.setNMonomer(nMonomer);

      // Allocate work space field arrays
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsRGrid_[i].allocate(dimensions);
         tmpFieldsKGrid_[i].allocate(dimensions);
      }
      isAllocatedRGrid_ = true;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocateFieldsBasis()
   {
      // Preconditions and constants
      UTIL_CHECK(hasMixture_);
      const int nMonomer = mixture_.nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = domain_.basis().nBasis();
      UTIL_CHECK(nBasis > 0);

      w_.allocateBasis(nBasis);
      c_.allocateBasis(nBasis);

      // Temporary work fields 
      tmpFieldsBasis_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsBasis_[i].allocate(nBasis);
      }
      isAllocatedBasis_ = true;
   }

   /*
   * Peek at field file header, initialize basis, allocate basis fields.
   */
   template <int D>
   void System<D>::allocateFieldsBasis(std::string filename)
   {
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(mixture_.nMonomer() > 0);

      // Open field file
      std::ifstream file;
      fileMaster_.openInputFile(filename, file);

      // Read field file header, and initialize basis if needed
      int nMonomer;
      UnitCell<D> unitCell;
      domain_.fieldIo().readFieldHeader(file, nMonomer, unitCell);
      // FieldIo::readFieldHeader initializes a basis if needed
      file.close();
      UTIL_CHECK(mixture_.nMonomer() == nMonomer);
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(domain_.basis().nBasis() > 0);

      // Allocate memory for fields
      allocateFieldsBasis(); 
   }

   /*
   * Read a filename string and echo to log file (used in readCommands).
   */
   template <int D>
   void System<D>::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      Log::file() << " " << Str(string, 20) << std::endl;
   }

   
   /*
   * Read and execute commands from a specified command file.
   */
   template <int D>
   void System<D>::readCommands(std::istream &in) 
   {
      UTIL_CHECK(isAllocatedRGrid_);
      std::string command, filename, inFileName, outFileName;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command <<std::endl;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W_BASIS") {
            readEcho(in, filename);
            readWBasis(filename);
         } else
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);
            readWRGrid(filename);
         } else
         if (command == "SET_UNIT_CELL") {
            UnitCell<D> unitCell;
            in >> unitCell;
            Log::file() << "   " << unitCell << std::endl;
            setUnitCell(unitCell);
         } else
         if (command == "READ_H_BASIS") {
            readEcho(in, filename);
            if (!h_.isAllocatedBasis()) {
               h_.allocateBasis(basis().nBasis());
            }
            if (!h_.isAllocatedRGrid()) {
               h_.allocateRGrid(mesh().dimensions());
            }
            h_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "READ_H_RGRID") {
            readEcho(in, filename);
            if (!h_.isAllocatedRGrid()) {
               h_.allocateRGrid(mesh().dimensions());
            }
            h_.readRGrid(filename, domain_.unitCell());
         } else
         if (command == "READ_MASK_BASIS") {
            UTIL_CHECK(domain_.basis().isInitialized());
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(basis().nBasis(), mesh().dimensions());
            }
            mask_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "READ_MASK_RGRID") {
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(basis().nBasis(), mesh().dimensions());
            }
            mask_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "COMPUTE") {
            // Solve the modified diffusion equation, without iteration
            compute();
         } else
         if (command == "ITERATE") {
            // Attempt iteration to convergence
            bool isContinuation = false;
            int fail = iterate(isContinuation);
            if (fail) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            // Do a series of iterations.
            sweep();
         } else
         if (command == "COMPRESS") {
            // Impose incompressibility
            UTIL_CHECK(hasCompressor());
            compressor().compress();
         } else
         if (command == "WRITE_W_BASIS") {
            readEcho(in, filename);
            writeWBasis(filename);
         } else 
         if (command == "WRITE_W_RGRID") {
            readEcho(in, filename);
            writeWRGrid(filename);
         } else 
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            writeCBasis(filename);
         } else
         if (command == "WRITE_C_RGRID") {
            readEcho(in, filename);
            writeCRGrid(filename);
         } else
         if (command == "WRITE_C_BLOCK_RGRID") {
            readEcho(in, filename);
            writeBlockCRGrid(filename);
         } else
         if (command == "WRITE_Q_SLICE") {
            int polymerId, blockId, directionId, segmentId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            in >> segmentId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n"
                        << Str("segment ID  ", 21) << segmentId << std::endl;
            writeQSlice(filename, polymerId, blockId, directionId, 
                                  segmentId);
         } else
         if (command == "WRITE_Q_TAIL") {
            readEcho(in, filename);
            int polymerId, blockId, directionId;
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQTail(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q") {
            readEcho(in, filename);
            int polymerId, blockId, directionId;
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(in, filename);
            writeQAll(filename);
         } else
         if (command == "WRITE_PARAM") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file);
            writeParamNoSweep(file);
            file.close();
         } else
         if (command == "WRITE_THERMO") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file, 
                                        std::ios_base::app);
            writeThermo(file);
            file.close();
         } else
         if (command == "WRITE_STARS") {
            readEcho(in, filename);
            writeStars(filename);
         } else
         if (command == "WRITE_WAVES") {
            readEcho(in, filename);
            writeWaves(filename);
         } else 
         if (command == "WRITE_MASK_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            UTIL_CHECK(mask_.isSymmetric());
            fieldIo().writeFieldBasis(filename, mask_.basis(), unitCell());
         } else 
         if (command == "WRITE_MASK_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            fieldIo().writeFieldRGrid(filename, mask_.rgrid(), unitCell());
         } else 
         if (command == "WRITE_H_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            UTIL_CHECK(h_.isSymmetric());
            fieldIo().writeFieldsBasis(filename, h_.basis(), unitCell());
         } else 
         if (command == "WRITE_H_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            fieldIo().writeFieldsRGrid(filename, h_.rgrid(), unitCell());
         } else 
         if (command == "BASIS_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            basisToRGrid(inFileName, outFileName);
         } else 
         if (command == "RGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            rGridToBasis(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToRGrid(inFileName, outFileName);
         } else
         if (command == "RGRID_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            rGridToKGrid(inFileName, outFileName);
         } else
         if (command == "BASIS_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            basisToKGrid(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToBasis(inFileName, outFileName);
         } else
         if (command == "CHECK_RGRID_SYMMETRY") {
            readEcho(in, inFileName);
            bool hasSymmetry;
            hasSymmetry = checkRGridFieldSymmetry(inFileName);
            if (hasSymmetry) {
               Log::file() 
                   << "Symmetry of r-grid file matches this space group." 
                   << std::endl;
            } else {
               Log::file() 
                   << "Symmetry of r-grid file does not match this space group" 
                   << std::endl
                   << "to within our error threshold of 1E-8."
                   << std::endl;
            }
         } else
         if (command == "GUESS_W_FROM_C") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            guessWfromC(inFileName, outFileName);
         } else
         if (command == "COMPARE_BASIS") {

            // Get two filenames for comparison
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            
            DArray< DArray<double> > Bfield1, Bfield2;
            domain_.fieldIo().readFieldsBasis(filecompare1, Bfield1, 
                                      domain_.unitCell());
            domain_.fieldIo().readFieldsBasis(filecompare2, Bfield2, 
                                      domain_.unitCell());
            // Note: Bfield1 & Bfield2 are allocated by readFieldsBasis

            // Compare and output report
            compare(Bfield1, Bfield2);

         } else
         if (command == "COMPARE_RGRID") {
            // Get two filenames for comparison
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            
            DArray< RField<D> > Rfield1, Rfield2;
            domain_.fieldIo().readFieldsRGrid(filecompare1, Rfield1, 
                                              domain_.unitCell());
            domain_.fieldIo().readFieldsRGrid(filecompare2, Rfield2, 
                                              domain_.unitCell());
            // Note: Rfield1 and Rfield2 will be allocated by readFieldsRGrid

            // Compare and output report
            compare(Rfield1, Rfield2);

         } else
         {
            Log::file() << "Error: Unknown command  " 
                        << command << std::endl;
            readNext = false;
         }
      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   template <int D>
   void System<D>::readCommands()
   {  
      if (fileMaster_.commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster_.commandFile()); 
   }

   // Chemical Potential Field Modifier Functions

   /*
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(filename); 
      }

      // Read w fields
      w_.readBasis(filename, domain_.unitCell());
      mixture_.setupUnitCell(domain_.unitCell());
      hasCFields_ = false;
   }

   /*
   * Read w-fields in real-space grid (r-grid) format.
   */
   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {

      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(filename); 
      }

      // Read w fields
      w_.readRGrid(filename, domain_.unitCell());
      mixture_.setupUnitCell(domain_.unitCell());
      hasCFields_ = false;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedBasis_);
      w_.setBasis(fields);
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray<Field> const & fields)
   {
      UTIL_CHECK(isAllocatedRGrid_);
      w_.setRGrid(fields);
      hasCFields_ = false;
   }

   // Unit Cell Modifier / Setter 

   /*
   * Set parameters of the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      domain_.setUnitCell(unitCell);
      mixture_.setupUnitCell(domain_.unitCell());
   }

   /*
   * Set state of the system unit cell.
   */
   template <int D>
   void 
   System<D>::setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                          FSArray<double, 6> const & parameters)
   {
      domain_.setUnitCell(lattice, parameters);
      mixture_.setupUnitCell(domain_.unitCell());
   }

   /*
   * Set parameters of the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      domain_.setUnitCell(parameters);
      mixture_.setupUnitCell(domain_.unitCell());
   }

   // Primary Field Theory Computations

   /*
   * Solve MDE for current w-fields, without iteration.
   */
   template <int D>
   void System<D>::compute(bool needStress)
   {
      UTIL_CHECK(w_.isAllocatedRGrid());
      UTIL_CHECK(c_.isAllocatedRGrid());
      UTIL_CHECK(w_.hasData());

      // Solve the modified diffusion equation (without iteration)
      mixture_.compute(w_.rgrid(), c_.rgrid(), mask_.phiTot());
      hasCFields_ = true;

      // Compute stress if requested
      if (needStress) {
         mixture_.computeStress();
      }

      // If w fields are symmetric, compute basis componens for c-fields
      if (w_.isSymmetric()) {
         UTIL_CHECK(c_.isAllocatedBasis());
         domain_.fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis());
      }

   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   template <int D>
   int System<D>::iterate(bool isContinuation)
   {
      UTIL_CHECK(iteratorPtr_);
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      hasCFields_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator (return 0 for convergence, 1 for failure)
      int error = iterator().solve(isContinuation);
      hasCFields_ = true;

      // If converged, compute related properties
      if (!error) {   
         if (!iterator().isFlexible()) {
            mixture().computeStress();
         }
         computeFreeEnergy();
         writeThermo(Log::file());
      }
      return error;
   }

   /*
   * Perform sweep along a line in parameter space.
   */
   template <int D>
   void System<D>::sweep()
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      UTIL_CHECK(hasSweep());
      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
   }
  
   // Field Theoretic Monte-Carlo Simulation
 
   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void System<D>::simulate(int nStep)
   {

      mcMoveManager().setup();
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for (int iStep = 0; iStep < nStep; ++iStep) {

         #if 0
         // Call analyzers
         if (Analyzer::baseInterval != 0) {
            if (iStep % Analyzer::baseInterval == 0) {
               if (analyzerManager().size() > 0) {
                  system().positionSignal().notify();
                  analyzerManager().sample(iStep);
                  system().positionSignal().notify();
               }
            }
         }
         #endif

         // Choose and attempt an McMove
         mcMoveManager().chooseMove().move();

      }
      timer.stop();
      double time = timer.time();

      #if 0
      // Final analyzers
      assert(iStep == nStep);
      if (Analyzer::baseInterval > 0) {
         if (iStep % Analyzer::baseInterval == 0) {
            if (analyzerManager().size() != 0) {
               system().positionSignal().notify();
               analyzerManager().sample(iStep);
               system().positionSignal().notify();
            }
         }
      }

      // Output results of all analyzers to output files
      if (Analyzer::baseInterval > 0) {
         analyzerManager().output();
      }
      #endif

      // Output results of move statistics to files
      mcMoveManager().output();

      // Output time for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep         " << nStep << std::endl;
      Log::file() << "run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep  " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Print McMove acceptance statistics
      long attempt;
      long accept;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(32) << left <<  "Move Name"
           << setw(12) << right << "Attempted"
           << setw(12) << right << "Accepted"
           << setw(15) << right << "AcceptRate"
           << endl;
      int nMove = mcMoveManager().size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = mcMoveManager()[iMove].nAttempt();
         accept  = mcMoveManager()[iMove].nAccept();
         Log::file() << setw(32) << left
              << mcMoveManager().className(iMove)
              << setw(12) << right << attempt
              << setw(12) << accept
              << setw(15) << fixed << setprecision(6)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

   /*
   * Save the current r-grid w fields to temporary storage. 
   *
   * Used before attempting a Monte-Carlo move.
   */
   template <int D>
   void System<D>::saveWRGrid()
   {
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(w_.hasData());

      int nMonomer = mixture_.nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsRGrid_[i] = w_.rgrid(i);
      }
   }

   /*
   * Restore save r-grid fields. 
   *
   * Used when an attempted Monte-Carlo move is rejected.
   */
   template <int D>
   void System<D>::restoreWRGrid()
   {  setWRGrid(tmpFieldsRGrid_); }

   // Thermodynamic Properties

   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());

      // Initialize to zero
      fHelmholtz_ = 0.0;
      fIdeal_ = 0.0;
      fInter_ = 0.0;
 
      double phi, mu;
      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      // Compute polymer ideal gas contributions to fHelhmoltz_
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 )/length;
            }
         }
      }

      // Compute solvent ideal gas contributions to fHelhmoltz_
      if (ns > 0) {
         Solvent<D>* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture_.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      int nm  = mixture_.nMonomer();
      int nBasis = domain_.basis().nBasis();

      double temp(0.0);
      // Compute Legendre transform subtraction
      // Use expansion in symmetry-adapted orthonormal basis
      for (int i = 0; i < nm; ++i) {
         for (int k = 0; k < nBasis; ++k) {
            temp -= w_.basis(i)[k] * c_.basis(i)[k];
         }
      }

      // If the system has a mask, then the volume that should be used
      // in calculating free energy/pressure is the volume available to
      // the polymers, not the total unit cell volume. We thus divide
      // all terms that involve integrating over the unit cell volume by
      // mask().phiTot(), the volume fraction of the unit cell that is 
      // occupied by the polymers. This properly scales them to the 
      // correct value. fExt_, fInter_, and the Legendre transform 
      // component of fIdeal_ all require this scaling. If no mask is 
      // present, mask.phiTot() = 1 and no scaling occurs.
      temp /= mask().phiTot(); 
      fIdeal_ += temp;
      fHelmholtz_ += fIdeal_;

      // Compute contribution from external fields, if fields exist
      if (hasExternalFields()) {
         fExt_ = 0.0;
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               fExt_ += h_.basis(i)[k] * c_.basis(i)[k];
            }
         }
         fExt_ /= mask().phiTot();
         fHelmholtz_ += fExt_;
      }

      // Compute excess interaction free energy [ phi^{T}*chi*phi ]
      double chi;
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
            chi = interaction().chi(i,j);
            for (int k = 0; k < nBasis; ++k) {
               fInter_ += chi * c_.basis(i)[k] * c_.basis(j)[k];
            }
         }
      }
      fInter_ /= mask().phiTot();
      fHelmholtz_ += fInter_;

      // Initialize pressure
      pressure_ = -fHelmholtz_;

      // Polymer corrections to pressure
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            if (phi > 1E-08) {
               pressure_ += mu * phi /length;
            }
         }
      }

      // Solvent corrections to pressure
      if (ns > 0) {
         Solvent<D>* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture_.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1E-08) {
               pressure_ += mu * phi /size;
            }
         }
      }

   }

   /*
   * Write thermodynamic properties to file.
   */
   template <int D>
   void System<D>::writeThermo(std::ostream& out) const
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      
      out << "Free energy components:" << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
      if (hasExternalFields()) {
         out << "fExt          " << Dbl(fExt_, 18, 11) << std::endl;
      }
      out << std::endl;

      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      if (np > 0) {
         out << "Polymers:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < np; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture_.polymer(i).phi(),18, 11)
                << "  " << Dbl(mixture_.polymer(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

      if (ns > 0) {
         out << "Solvents:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < ns; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture_.solvent(i).phi(),18, 11)
                << "  " << Dbl(mixture_.solvent(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

   }

   /*
   * Initialize the homogeneous_ member object, which describes 
   * thermodynamics of a homogeneous reference system.
   */
   template <int D>
   void System<D>::initHomogeneous()
   {
      UTIL_CHECK(hasMixture_);

      // Set number of molecular species and monomers
      int nm = mixture_.nMonomer(); 
      int np = mixture_.nPolymer(); 
      int ns = mixture_.nSolvent(); 
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      int i;   // molecule index
      int j;   // monomer index
 
      // Loop over polymer molecule species
      if (np > 0) {

         // Allocate array of clump sizes
         DArray<double> s;
         s.allocate(nm);
      
         int k;   // block or clump index
         int nb;  // number of blocks
         int nc;  // number of clumps

         // Loop over polymer species
         for (i = 0; i < np; ++i) {
   
            // Initial array of clump sizes for this polymer
            for (j = 0; j < nm; ++j) {
               s[j] = 0.0;
            }
   
            // Compute clump sizes for all monomer types.
            nb = mixture_.polymer(i).nBlock(); 
            for (k = 0; k < nb; ++k) {
               Block<D>& block = mixture_.polymer(i).block(k);
               j = block.monomerId();
               s[j] += block.length();
            }
    
            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (s[j] > 1.0E-8) {
                  ++nc;
               }
            }
            homogeneous_.molecule(i).setNClump(nc);
    
            // Set clump properties for this Homogeneous::Molecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (s[j] > 1.0E-8) {
                  homogeneous_.molecule(i).clump(k).setMonomerId(j);
                  homogeneous_.molecule(i).clump(k).setSize(s[j]);
                  ++k;
               }
            }
            homogeneous_.molecule(i).computeSize();
   
         }
      }

      // Add solvent contributions
      if (ns > 0) {
         double size;
         int monomerId;
         for (int is = 0; is < ns; ++is) {
            i = is + np;
            monomerId = mixture_.solvent(is).monomerId();
            size = mixture_.solvent(is).size();
            homogeneous_.molecule(i).setNClump(1);
            homogeneous_.molecule(i).clump(0).setMonomerId(monomerId);
            homogeneous_.molecule(i).clump(0).setSize(size);
            homogeneous_.molecule(i).computeSize();
         }
      }

   }

   // Output functions

   /*
   * Write w-fields in symmetry-adapted basis format. 
   */
   template <int D>
   void System<D>::writeWBasis(std::string const & filename) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      domain_.fieldIo().writeFieldsBasis(filename, w_.basis(), 
                                         domain_.unitCell());
   }

   /*
   * Write w-fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(w_.hasData());
      domain_.fieldIo().writeFieldsRGrid(filename, w_.rgrid(), 
                                         domain_.unitCell());
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(const std::string & filename) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(w_.isSymmetric());
      domain_.fieldIo().writeFieldsBasis(filename, c_.basis(), 
                                         domain_.unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(hasCFields_);
      domain_.fieldIo().writeFieldsRGrid(filename, c_.rgrid(), 
                                         domain_.unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format, for 
   * each block (or solvent) individually rather than for each species.
   */
   template <int D>
   void System<D>::writeBlockCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedRGrid_);
      UTIL_CHECK(hasCFields_);

      // Create and allocate the DArray of fields to be written
      DArray<Field> blockCFields;
      blockCFields.allocate(mixture_.nSolvent() + mixture_.nBlock());
      int n = blockCFields.capacity();
      for (int i = 0; i < n; i++) {
         blockCFields[i].allocate(domain_.mesh().dimensions());
      }

      // Get data from Mixture and write to file
      mixture_.createBlockCRGrid(blockCFields);
      domain_.fieldIo().writeFieldsRGrid(filename, blockCFields, 
                                         domain_.unitCell());
   }

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQSlice(const std::string & filename, 
                               int polymerId, int blockId, 
                               int directionId, int segmentId) 
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator<D> const & 
           propagator = polymer.propagator(blockId, directionId);
      RField<D> const& field = propagator.q(segmentId);
      domain_.fieldIo().writeFieldRGrid(filename, field, 
                                        domain_.unitCell());
   }

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQTail(const std::string & filename, 
                              int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      RField<D> const& 
            field = polymer.propagator(blockId, directionId).tail();
      domain_.fieldIo().writeFieldRGrid(filename, field, 
                                        domain_.unitCell());
   }

   /*
   * Write the propagator for a block and direction.
   */
   template <int D>
   void System<D>::writeQ(const std::string & filename, 
                          int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator<D> const& propagator 
                              = polymer.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);

      // Write header
      fieldIo().writeFieldHeader(file, 1, domain_.unitCell());
      file << "ngrid" << std::endl
           << "          " << domain_.mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
          file << "slice " << i << std::endl;
          fieldIo().writeFieldRGrid(file, propagator.q(i), 
                                    domain_.unitCell(), hasHeader);
      }
   }

   /*
   * Write propagators for all blocks of all polymers to files.
   */
   template <int D>
   void System<D>::writeQAll(std::string const & basename)
   {
      std::string filename;
      int np, nb, ip, ib, id;
      np = mixture_.nPolymer();
      for (ip = 0; ip < np; ++ip) {
         //Polymer<D> const * polymerPtr = &mixture_.polymer(ip);
         //nb = polymerPtr->nBlock();
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
               filename += ".rf";
               writeQ(filename, ip, ib, id);
            }
         }
      }
   }

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D>
   void System<D>::writeStars(const std::string & outFileName) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(),
                                 domain_.unitCell());
      domain_.basis().outputStars(outFile);
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D>
   void System<D>::writeWaves(const std::string & outFileName) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(), 
                                 domain_.unitCell());
      domain_.basis().outputWaves(outFile);
   }

   // Field conversion command functions

   /*
   * Convert fields from symmetry-adpated basis to real-space grid format.
   */
   template <int D>
   void System<D>::basisToRGrid(const std::string & inFileName,
                                const std::string & outFileName)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert, and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToRGrid(tmpFieldsBasis_, tmpFieldsRGrid_);
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_, 
                                 tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::rGridToBasis(const std::string & inFileName,
                                const std::string & outFileName) 
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      fieldIo().convertRGridToBasis(tmpFieldsRGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
   */
   template <int D>
   void System<D>::kGridToRGrid(const std::string & inFileName,
                                const std::string& outFileName)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         domain_.fft().inverseTransform(tmpFieldsKGrid_[i], 
                                        tmpFieldsRGrid_[i]);
      }
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_, 
                                 tmpUnitCell);
   }

   /*
   * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
   */
   template <int D>
   void System<D>::rGridToKGrid(const std::string & inFileName,
                                const std::string & outFileName)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, 
                                tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         domain_.fft().forwardTransform(tmpFieldsRGrid_[i], 
                                        tmpFieldsKGrid_[i]);
      }
      domain_.fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, 
                                         tmpUnitCell);
   }

   /*
   * Convert fields from Fourier (k-grid) to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::kGridToBasis(const std::string & inFileName,
                                const std::string& outFileName)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      domain_.fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, 
                                        tmpUnitCell);
      domain_.fieldIo().convertKGridToBasis(tmpFieldsKGrid_, 
                                            tmpFieldsBasis_);
      domain_.fieldIo().writeFieldsBasis(outFileName, 
                                         tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert fields from symmetry-adapted basis to Fourier (k-grid) format.
   */
   template <int D>
   void System<D>::basisToKGrid(const std::string & inFileName,
                                const std::string & outFileName) 
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      domain_.fieldIo().readFieldsBasis(inFileName, 
                                        tmpFieldsBasis_, tmpUnitCell);
      domain_.fieldIo().convertBasisToKGrid(tmpFieldsBasis_, 
                                            tmpFieldsKGrid_);
      domain_.fieldIo().writeFieldsKGrid(outFileName,  
                                         tmpFieldsKGrid_, tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   bool System<D>::checkRGridFieldSymmetry(const std::string & inFileName) 
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }

      // Read fields
      UnitCell<D> tmpUnitCell;
      domain_.fieldIo().readFieldsRGrid(inFileName, 
                                        tmpFieldsRGrid_, tmpUnitCell);

      // Check symmetry for all fields
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         bool symmetric;
         symmetric = domain_.fieldIo().hasSymmetry(tmpFieldsRGrid_[i]);
         if (!symmetric) {
            return false;
         }
      }
      return true;

   }

   /*
   * Compare two fields in basis format.
   */ 
   template <int D>
   void System<D>::compare(const DArray< DArray<double> > field1, 
                           const DArray< DArray<double> > field2)
   {
      BFieldComparison comparison(1);
      comparison.compare(field1,field2);

      Log::file() << "\n Basis expansion field comparison results" 
                  << std::endl;
      Log::file() << "     Maximum Absolute Difference:   " 
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   " 
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   /*
   * Compare two fields in coordinate grid format.
   */ 
   template <int D>
   void System<D>::compare(const DArray< RField<D> > field1, 
                           const DArray< RField<D> > field2)
   {
      RFieldComparison<D> comparison;
      comparison.compare(field1, field2);

      Log::file() << "\n Real-space field comparison results" 
                  << std::endl;
      Log::file() << "     Maximum Absolute Difference:   " 
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   " 
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   /*
   * Construct guess for omega (w-field) from rho (c-field).
   *
   * Modifies wFields and wFieldsRGrid and outputs wFields.
   */
   template <int D>
   void System<D>::guessWfromC(std::string const & inFileName, 
                               std::string const & outFileName)
   {
      const int nm = mixture_.nMonomer();
      UTIL_CHECK(nm > 0);

      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         allocateFieldsBasis(inFileName); 
      }
      const int nb = domain_.basis().nBasis();
      UTIL_CHECK(nb > 0);

      // Read c fields and set unit cell
      domain_.fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, 
                                        domain_.unitCell());

      DArray<double> wtmp;
      wtmp.allocate(nm);

      // Compute w fields from c fields
      int i, j, k;
      for (i = 0; i < nb; ++i) {
         for (j = 0; j < nm;  ++j) {
            wtmp[j] = 0.0;
            for (k = 0; k < nm; ++k) {
               wtmp[j] += interaction().chi(j,k)*tmpFieldsBasis_[k][i];
            }
         }
         for (j = 0; j < nm;  ++j) {
            tmpFieldsBasis_[j][i] = wtmp[j];
         }
      }

      // Set initial guess, and write to file
      w_.setBasis(tmpFieldsBasis_);
      domain_.fieldIo().writeFieldsBasis(outFileName, w_.basis(), 
                                         domain_.unitCell());

      hasCFields_ = false;
   }

} // namespace Pspc
} // namespace Pscf
#endif
