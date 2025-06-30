#ifndef RPC_SYSTEM_TPP
#define RPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/simulator/SimulatorFactory.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/scft/sweep/Sweep.h>
#include <rpc/scft/sweep/SweepFactory.h>
#include <rpc/scft/iterator/Iterator.h>
#include <rpc/scft/iterator/IteratorFactory.h>
#include <rpc/environment/EnvironmentFactory.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>


#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>

#include <pscf/environment/Environment.h>
#include <pscf/inter/Interaction.h>
#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/param/BracketPolicy.h>
#include <util/param/ParamComponent.h>
#include <util/signal/Signal.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

#include <string>
#include <unistd.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   System<D>::System()
    : mixture_(),
      domain_(),
      fileMaster_(),
      w_(),
      c_(),
      h_(),
      mask_(),
      interactionPtr_(nullptr),
      environmentPtr_(nullptr),
      environmentFactoryPtr_(nullptr),
      iteratorPtr_(nullptr),
      iteratorFactoryPtr_(nullptr),
      sweepPtr_(nullptr),
      sweepFactoryPtr_(nullptr),
      simulatorPtr_(nullptr),
      simulatorFactoryPtr_(nullptr),
      fHelmholtz_(0.0),
      fIdeal_(0.0),
      fInter_(0.0),
      fExt_(0.0),
      pressure_(0.0),
      polymerModel_(PolymerModel::Thread),
      isAllocatedGrid_(false),
      isAllocatedBasis_(false),
      hasMixture_(false),
      hasCFields_(false),
      hasFreeEnergy_(false),
      hasStress_(false)
   {
      setClassName("System");  // Set block label in parameter file
      BracketPolicy::set(BracketPolicy::Optional);


      // Create associations among class members
      domain_.setFileMaster(fileMaster_);
      w_.setFieldIo(domain_.fieldIo());
      w_.setReadUnitCell(domain_.unitCell());  
      w_.setWriteUnitCell(domain_.unitCell());
      h_.setFieldIo(domain_.fieldIo());        
      h_.setReadUnitCell(tmpUnitCell_);       
      h_.setWriteUnitCell(domain_.unitCell());
      mask_.setFieldIo(domain_.fieldIo());
      mask_.setReadUnitCell(tmpUnitCell_);       
      mask_.setWriteUnitCell(domain_.unitCell());
      c_.setFieldIo(domain_.fieldIo());
      c_.setWriteUnitCell(domain_.unitCell());

      // Note: When w_ is read from a file  in basis or r-grid format,
      // the parameters of the system unit cell, domain_.unitCell(), are
      // set to those in the field file header. When imposed fields h_
      // and mask_ are read from a file, however unit cell parameters 
      // from the field file header are read into a mutable workspace 
      // object, tmpUnitCell_, and ignored.

      // Create dynamically allocated objects owned by this System
      interactionPtr_ = new Interaction();
      environmentFactoryPtr_ = new EnvironmentFactory<D>(*this);
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this);
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      simulatorFactoryPtr_ = new SimulatorFactory<D>(*this);

      // Use of "signals" to maintain data relationships:
      // Signals are instances of Unit::Signal<void> that are used to 
      // notify "observer" objects of modification of data owned by a 
      // related "notifier" object. Each signal is owned by a notifier
      // object that maintains data that may be modified. Each signal
      // maintains a list of observers objects that should be notified
      // whenever the data owned by the notifier object changes. Each 
      // observer is added by the Signal<void>::addObserver function
      // template, which takes two arguments: a reference to an observer 
      // object (i.e., an instance of a class) and a pointer to a member 
      // function of that class which will be invoked on that object when 
      // the signal is triggered by modification of associated data.

      // Addition of observers to signals

      // Signal triggered by unit cell modification
      Signal<void>& cellSignal = domain_.unitCell().signal();
      cellSignal.addObserver(*this, &System<D>::clearUnitCellData);

      // Signal triggered by basis construction
      Signal<void>& basisSignal = domain_.basis().signal();
      basisSignal.addObserver(*this, &System<D>::allocateFieldsBasis);

      // Signal triggered by w-field modification
      w_.signal().addObserver(*this, &System<D>::clearCFields);

      // Signal triggered by h-field modification
      h_.signal().addObserver(*this, &System<D>::clearCFields);

      // Signal triggered by mask modification
      mask_.signal().addObserver(*this, &System<D>::clearCFields);
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
      if (environmentPtr_) {
         delete environmentPtr_;
      }
      if (environmentFactoryPtr_) {
         delete environmentFactoryPtr_;
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
      if (simulatorPtr_) {
         delete simulatorPtr_;
      }
      if (simulatorFactoryPtr_) {
         delete simulatorFactoryPtr_;
      }
   }

   // Lifetime (called in main program)

   /*
   * Process command line options.
   */
   template <int D>
   void System<D>::setOptions(int argc, char **argv)
   {
      bool eFlag = false;  // echo
      bool dFlag = false;  // spatial dimension (1, 2, or 3)
      bool pFlag = false;  // param file
      bool cFlag = false;  // command file
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      bool tFlag = false;  // nThread
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
      int dArg = 0;
      int tArg = 0;

      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "ed:p:c:i:o:t:")) != -1) {
         switch (c) {
         case 'e':
            eFlag = true;
            break;
         case 'd':
            dFlag = true;
            dArg = atoi(optarg);
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
         case 't': // number of threads, user set
            tFlag = true;
            tArg = atoi(optarg);
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }

      // Set flag to echo parameters as they are read.
      if (eFlag) {
         Util::ParamComponent::setEcho(true);
      }

      // Check -d flag
      if (dFlag) {
         UTIL_CHECK(D == dArg);
      } else {
         UTIL_THROW("Missing required -d option");
      }

      // Check option -p, set parameter file name
      if (pFlag) {
         fileMaster_.setParamFileName(std::string(pArg));
      } else {
         UTIL_THROW("Missing required -p option - no parameter file");
      }

      // Check option -c, set command file name
      if (cFlag) {
         fileMaster_.setCommandFileName(std::string(cArg));
      } else {
         UTIL_THROW("Missing required -c option - no command file");
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster_.setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster_.setOutputPrefix(std::string(oArg));
      }

      if (tFlag) {
         if (tArg <= 0) {
            UTIL_THROW("Error: Non-positive thread count -t option");
         }
      }
   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {
      // Optionally read polymerModel_ enum value, set the global value
      if (!PolymerModel::isLocked()) {
         polymerModel_ = PolymerModel::Thread;
         readOptional(in, "polymerModel", polymerModel_);

         // Set the global enumeration value
         PolymerModel::setModel(polymerModel_);
      }

      // Check number of times the global polymer model has been set.
      // Retain this value so we can check at the end of this function
      // that it was not set again within the function.

      int nSetPolymerModel = PolymerModel::nSet();

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

      // Read the Interaction{ ... } block
      interaction().setNMonomer(nm);
      readParamComposite(in, interaction());

      // Read the Domain{ ... } block
      readParamComposite(in, domain_);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(domain_.unitCell().nParameter() > 0);
      UTIL_CHECK(domain_.unitCell().lattice() != UnitCell<D>::Null);
      domain_.fieldIo().setNMonomer(nm);

      // Setup the mixture
      mixture_.associate(domain_.mesh(), domain_.fft(),
                         domain_.unitCell(), domain_.waveList());
      mixture_.allocate();

      // Allocate memory for w and c fields in r-grid form
      allocateFieldsGrid();

      // Optionally construct an Environment object
      std::string className;
      bool isEnd;
      environmentPtr_ =
         environmentFactoryPtr_->readObjectOptional(in, *this, className,
                                                    isEnd);
      if (!environmentPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Environment{ [absent] }\n";
      }

      // Optionally construct an Iterator object
      if (!isEnd) {
         iteratorPtr_ =
            iteratorFactoryPtr_->readObjectOptional(in, *this, className,
                                                   isEnd);
         if (!iteratorPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Iterator{ [absent] }\n";
         }
      }

      // Optionally construct a Sweep object (if an iterator exists)
      if (hasIterator() && !isEnd) {
         sweepPtr_ =
            sweepFactoryPtr_->readObjectOptional(in, *this, className,
                                                 isEnd);
         if (!sweepPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Sweep{ [absent] }\n";
         }
      }

      // Optionally construct a Simulator object
      if (!isEnd) {
         simulatorPtr_ =
            simulatorFactoryPtr_->readObjectOptional(in, *this,
                                                     className, isEnd);
         if (!simulatorPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Simulator{ [absent] }\n";
         }
      }

      // Check that the polymer model was not reset after initialization
      UTIL_CHECK(PolymerModel::nSet() == nSetPolymerModel);
   }

   /*
   * Read parameter file (including open and closing brackets).
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
   * Read and execute commands from a specified command file.
   */
   template <int D>
   void System<D>::readCommands(std::istream &in)
   {
      UTIL_CHECK(isAllocatedGrid_);
      std::string command, filename, inFileName, outFileName;
      FieldIo<D> const & fieldIo = domain_.fieldIo();

      bool readNext = true;
      while (readNext) {

         in >> command;

         if (in.eof()) {
            break;
         } else {
            Log::file() << command << std::endl;
         }

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W_BASIS") {
            readEcho(in, filename);
            w_.readBasis(filename);
            UTIL_CHECK(domain_.unitCell().isInitialized());
            UTIL_CHECK(domain_.basis().isInitialized());
            UTIL_CHECK(!domain_.waveList().hasKSq());
            UTIL_CHECK(isAllocatedBasis_);
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(c_.hasData() == hasCFields_);
            UTIL_CHECK(!hasFreeEnergy_);
            UTIL_CHECK(!hasStress_);
         } else
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);
            w_.readRGrid(filename);
            UTIL_CHECK(domain_.unitCell().isInitialized());
            UTIL_CHECK(!domain_.waveList().hasKSq());
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(!hasFreeEnergy_);
            UTIL_CHECK(!hasStress_);
            UTIL_CHECK(c_.hasData() == hasCFields_);
         } else
         if (command == "ESTIMATE_W_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.estimateWBasis(inFileName, outFileName, 
                                   interaction().chi());
         } else
         if (command == "SET_UNIT_CELL") {
            UnitCell<D> unitCell;
            in >> unitCell;
            Log::file() << "   " << unitCell << std::endl;
            setUnitCell(unitCell);
         } else
         if (command == "COMPUTE") {
            // Solve the modified diffusion equation, without iteration
            compute();
         } else
         if (command == "ITERATE") {
            // Attempt to iteratively solve a single SCFT problem
            bool isContinuation = false;
            int fail = iterate(isContinuation);
            if (fail) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            // Attempt to solve a sequence of SCFT problems along a path
            // through parameter space
            sweep();
         } else
         if (command == "COMPRESS") {
            // Impose incompressibility
            UTIL_CHECK(hasSimulator());
            simulator().compressor().compress();
         } else
         if (command == "SIMULATE") {
            // Perform a field theoretic simulation
            int nStep;
            in >> nStep;
            Log::file() << "   "  << nStep << "\n";
            simulate(nStep);
         } else
         if (command == "ANALYZE" || command == "ANALYZE_TRAJECTORY") {
            // Read and analyze a field trajectory file
            int min, max;
            in >> min;
            in >> max;
            Log::file() << "   "  << min ;
            Log::file() << "   "  << max << "\n";
            std::string classname;
            readEcho(in, classname);
            readEcho(in, filename);
            simulator().analyze(min, max, classname, filename);
         } else
         if (command == "WRITE_PARAM") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster_.openOutputFile(filename, file);
            writeParamNoSweep(file);
            file.close();
         } else
         if (command == "WRITE_THERMO") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster_.openOutputFile(filename, file,
                                       std::ios_base::app);
            writeThermo(file);
            file.close();
         } else
         if (command == "WRITE_STRESS") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster_.openOutputFile(filename, file,
                                       std::ios_base::app);
            writeStress(file);
            file.close();
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
         if (command == "WRITE_BLOCK_C_RGRID") {
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
                        << Str("direction ID  ", 21) << directionId
                        << "\n"
                        << Str("segment ID  ", 21) << segmentId
                        << std::endl;
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
                        << Str("direction ID  ", 21) << directionId
                        << "\n";
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
                        << Str("direction ID  ", 21) << directionId
                        << "\n";
            writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(in, filename);
            writeQAll(filename);
         } else
         if (command == "WRITE_STARS") {
            readEcho(in, filename);
            domain_.writeStars(filename);
         } else
         if (command == "WRITE_WAVES") {
            readEcho(in, filename);
            domain_.writeWaves(filename);
         } else
         if (command == "WRITE_GROUP") {
            readEcho(in, filename);
            domain_.writeGroup(filename);
         } else
         if (command == "BASIS_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertBasisToRGrid(inFileName, outFileName);
         } else
         if (command == "RGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertRGridToBasis(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertKGridToRGrid(inFileName, outFileName);
         } else
         if (command == "RGRID_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertRGridToKGrid(inFileName, outFileName);
         } else
         if (command == "BASIS_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertBasisToKGrid(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.convertKGridToBasis(inFileName, outFileName);
         } else
         if (command == "CHECK_RGRID_SYMMETRY") {
            double epsilon;
            readEcho(in, inFileName);
            readEcho(in, epsilon);
            bool hasSymmetry;
            hasSymmetry = fieldIo.hasSymmetry(inFileName, epsilon);
            if (hasSymmetry) {
               Log::file() << std::endl
                   << "Symmetry of r-grid file matches this space group"
                   << std::endl << std::endl;
            } else {
               Log::file() << std::endl
                   << "Symmetry of r-grid file does not match this\n"
                   << "space group to within error threshold of "
                   << Dbl(epsilon) << " ." << std::endl << std::endl;
            }
         } else
         if (command == "COMPARE_BASIS") {
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            fieldIo.compareFieldsBasis(filecompare1, filecompare2);
         } else
         if (command == "COMPARE_RGRID") {
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            fieldIo.compareFieldsRGrid(filecompare1, filecompare2);
         } else
         if (command == "SCALE_BASIS") {
            double factor;
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            readEcho(in, factor);
            fieldIo.scaleFieldsBasis(inFileName, outFileName, factor);
         } else
         if (command == "SCALE_RGRID") {
            double factor;
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            readEcho(in, factor);
            fieldIo.scaleFieldsRGrid(inFileName, outFileName, factor);
         } else
         if (command == "EXPAND_RGRID_DIMENSION") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);

            // Read expanded spatial dimension d
            int d;
            in >> d;
            UTIL_CHECK(d > D);
            Log::file() << Str("Expand dimension to: ") << d << "\n";

            // Read numbers of grid points along added dimensions
            DArray<int> newGridDimensions;
            newGridDimensions.allocate(d-D);
            for (int i = 0; i < d-D; i++) {
               in >> newGridDimensions[i];
            }
            Log::file()
               << Str("Numbers of grid points in added dimensions :  ");
            for (int i = 0; i < d-D; i++) {
               Log::file() << newGridDimensions[i] << " ";
            }
            Log::file() << "\n";

            fieldIo.expandRGridDimension(inFileName, outFileName,
                                         d, newGridDimensions);

         } else
         if (command == "REPLICATE_UNIT_CELL") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);

            // Read numbers of replicas along each direction
            IntVec<D> replicas;
            for (int i = 0; i < D; i++){
               in >> replicas[i];
            }
            for (int i = 0; i < D; i++){
               Log::file() << "Replicate unit cell in direction "
                           << i << " : ";
               Log::file() << replicas[i] << " times ";
               Log::file() << "\n";
            }
            fieldIo.replicateUnitCell(inFileName, outFileName, replicas);

         } else
         if (command == "READ_H_BASIS") {
            readEcho(in, filename);
            h_.readBasis(filename);
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(c_.hasData() == hasCFields_);
         } else
         if (command == "READ_H_RGRID") {
            readEcho(in, filename);
            h_.readRGrid(filename);
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(c_.hasData() == hasCFields_);
         } else
         if (command == "WRITE_H_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            UTIL_CHECK(h_.isSymmetric());
            fieldIo.writeFieldsBasis(filename, h_.basis(),
                                     domain_.unitCell());
         } else
         if (command == "WRITE_H_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            fieldIo.writeFieldsRGrid(filename, h_.rgrid(),
                                     domain_.unitCell());
         } else
         if (command == "READ_MASK_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(domain_.basis().isInitialized());
            if (!mask_.isAllocatedBasis()) {
               mask_.allocateBasis(domain_.basis().nBasis());
            }
            if (!mask_.isAllocatedRGrid()) {
               mask_.allocateRGrid(domain_.mesh().dimensions());
            }
            mask_.readBasis(filename);
         } else
         if (command == "READ_MASK_RGRID") {
            readEcho(in, filename);
            if (!mask_.isAllocatedRGrid()) {
               mask_.allocateRGrid(domain_.mesh().dimensions());
            }
            if (iterator().isSymmetric() && !mask_.isAllocatedBasis()) {
               UTIL_CHECK(domain_.basis().isInitialized());
               mask_.allocateBasis(domain_.basis().nBasis());
            }
            mask_.readRGrid(filename);
         } else
         if (command == "WRITE_MASK_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            UTIL_CHECK(mask_.isSymmetric());
            fieldIo.writeFieldBasis(filename, mask_.basis(),
                                    domain_.unitCell());
         } else
         if (command == "WRITE_MASK_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            fieldIo.writeFieldRGrid(filename, mask_.rgrid(),
                                    domain_.unitCell(),
                                    mask_.isSymmetric());
         } else
         if (command == "WRITE_TIMERS") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster_.openOutputFile(filename, file);
            writeTimers(file);
            file.close();
         } else
         if (command == "CLEAR_TIMERS") {
            clearTimers();
         } else {
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

   // Unit Cell Modifiers

   /*
   * Set the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      // Preconditions
      UTIL_CHECK(domain_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(domain_.lattice() == unitCell.lattice());

      // Set system unit cell
      domain_.unitCell() = unitCell;

      // If necessary, make basis
      if (domain_.hasGroup() && !domain_.basis().isInitialized()) {
         domain_.makeBasis();
      }

      // Postconditions
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(domain_.unitCell().lattice() == domain_.lattice());
      if (domain_.hasGroup()) {
         UTIL_CHECK(domain_.basis().isInitialized());
         UTIL_CHECK(isAllocatedBasis_);
      }
      UTIL_CHECK(!domain_.waveList().hasKSq());
   }

   /*
   * Set parameters of the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      // Precondition
      UTIL_CHECK(domain_.lattice() != UnitCell<D>::Null);

      // Set system unit cell
      if (domain_.unitCell().lattice() == UnitCell<D>::Null) {
         domain_.unitCell().set(domain_.lattice(), parameters);
      } else {
         domain_.unitCell().setParameters(parameters);
      }

      // If necessary, make basis
      if (domain_.hasGroup() && !domain_.basis().isInitialized()) {
         domain_.makeBasis();
      }

      // Postconditions
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(domain_.unitCell().lattice() == domain_.lattice());
      if (domain_.hasGroup()) {
         UTIL_CHECK(domain_.basis().isInitialized());
         UTIL_CHECK(isAllocatedBasis_);
      }
      UTIL_CHECK(!domain_.waveList().hasKSq());
   }

   /*
   * Notify System members that unit cell parameters have been modified.
   */
   template <int D>
   void System<D>::clearUnitCellData()
   {
      clearCFields();
      mixture_.clearUnitCellData();
      domain_.waveList().clearUnitCellData();
      if (hasEnvironment()) {
         environment().reset();
      }
   }

   // Primary Field Theory Computations

   /*
   * Solve MDE for current w fields, without iteration.
   */
   template <int D>
   void System<D>::compute(bool needStress)
   {
      // Preconditions
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(w_.isAllocatedRGrid());
      UTIL_CHECK(c_.isAllocatedRGrid());
      UTIL_CHECK(w_.hasData());
      clearCFields();

      // Make sure Environment is up-to-date
      if (hasEnvironment()) {
         if (environment().needsUpdate()) {
            environment().generate();
         }
      }

      // Solve the modified diffusion equation (without iteration)
      mixture_.compute(w_.rgrid(), c_.rgrid(), mask_.phiTot());
      c_.setHasData(true);
      hasCFields_ = true;
      hasFreeEnergy_ = false;
      hasStress_ = false;

      // If w fields are symmetric, compute basis components for c fields
      if (w_.isSymmetric()) {
         UTIL_CHECK(c_.isAllocatedBasis());
         domain_.fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis(),
                                               false);
         c_.setIsSymmetric(true);
      }

      // Compute stress if needed
      if (needStress) {
         computeStress();
      }
   }

   /*
   * Iteratively solve a SCFT problem.
   */
   template <int D>
   int System<D>::iterate(bool isContinuation)
   {
      // Preconditions
      UTIL_CHECK(hasIterator());
      UTIL_CHECK(w_.hasData());
      if (iterator().isSymmetric()) {
         UTIL_CHECK(isAllocatedBasis_);
         UTIL_CHECK(w_.isSymmetric());
      }
      clearCFields();

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Make sure Environment is up-to-date
      if (hasEnvironment()) {
         if (environment().needsUpdate()) {
            environment().generate();
         }
      }

      // Call iterator (return 0 for convergence, 1 for failure)
      int error = iterator().solve(isContinuation);
      UTIL_CHECK(c_.hasData());
      UTIL_CHECK(c_.hasData() == hasCFields_);

      // If converged, compute related thermodynamic properties
      if (!error) {
         computeFreeEnergy(); // Sets hasFreeEnergy_ = true
         writeThermo(Log::file());
         if (!iterator().isFlexible()) {
            if (!mixture_.hasStress()) {
               mixture_.computeStress(mask().phiTot());
            }
            writeStress(Log::file());
         }
      }

      return error;
   }

   /*
   * Perform an SCFT sweep along a path in parameter space.
   */
   template <int D>
   void System<D>::sweep()
   {
      // Preconditions
      UTIL_CHECK(hasIterator());
      UTIL_CHECK(hasSweep());
      UTIL_CHECK(w_.hasData());
      if (iterator().isSymmetric()) {
         UTIL_CHECK(isAllocatedBasis_);
         UTIL_CHECK(w_.isSymmetric());
      }

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
   }

   /*
   * Perform a field theoretic simulation of nStep steps.
   */
   template <int D>
   void System<D>::simulate(int nStep)
   {
      UTIL_CHECK(nStep > 0);
      UTIL_CHECK(hasSimulator());
      clearCFields();

      simulator().simulate(nStep);
   }

   // SCFT Thermodynamic Properties

   /*
   * Compute Helmholtz free energy and pressure.
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      if (hasFreeEnergy_) return;

      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(c_.hasData());
      UTIL_CHECK(c_.hasData() == hasCFields_);

      if (hasEnvironment()) UTIL_CHECK(!environment().needsUpdate());

      int nm = mixture_.nMonomer();   // number of monomer types
      int np = mixture_.nPolymer();   // number of polymer species
      int ns = mixture_.nSolvent();   // number of solvent species

      // Initialize all free energy contributions to zero
      fHelmholtz_ = 0.0;
      fIdeal_ = 0.0;
      fInter_ = 0.0;
      fExt_ = 0.0;

      double phi, mu;

      // Compute polymer ideal gas contributions to fIdeal_
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            if (PolymerModel::isThread()) {
               length = polymerPtr->length();
            } else {
               length = (double) polymerPtr->nBead();
            }
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 ) / length;
               // Recall: mu = ln(phi/q)
            }
         }
      }

      // Compute solvent ideal gas contributions to fIdeal_
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

      // Volume integrals with a mask: If the system has a mask, then the
      // volume that should be used in calculating free energy/pressure
      // is the volume available to the material, not the total unit cell
      // volume. We thus divide all terms that involve integrating over
      // the unit cell volume by quantity mask().phiTot(), which is the
      // volume fraction of the unit cell that is occupied by material.
      // This properly scales them to the correct value. fExt_, fInter_,
      // and the Legendre transform component of fIdeal_ all require
      // this scaling. If no mask is present, mask.phiTot() = 1 and no
      // scaling occurs.
      const double phiTot = mask().phiTot();

      // Compute Legendre transform subtraction from fIdeal_
      double temp = 0.0;
      if (w_.isSymmetric()) {
         // Use expansion in symmetry-adapted orthonormal basis
         UTIL_CHECK(w_.isAllocatedBasis());
         UTIL_CHECK(c_.isAllocatedBasis());
         const int nBasis = domain_.basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               temp -= w_.basis(i)[k] * c_.basis(i)[k];
            }
         }
      } else {
         // Use summation over grid points
         const int meshSize = domain_.mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < meshSize; ++k) {
               temp -= w_.rgrid(i)[k] * c_.rgrid(i)[k];
            }
         }
         temp /= double(meshSize);
      }
      temp /= phiTot;
      fIdeal_ += temp;
      fHelmholtz_ += fIdeal_;

      // Compute contribution from external fields, if they exist
      if (hasExternalFields()) {
         if (w_.isSymmetric() && h_.isSymmetric()) {
            // Use expansion in symmetry-adapted orthonormal basis
            UTIL_CHECK(h_.isAllocatedBasis());
            UTIL_CHECK(c_.isAllocatedBasis());
            const int nBasis = domain_.basis().nBasis();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < nBasis; ++k) {
                  fExt_ += h_.basis(i)[k] * c_.basis(i)[k];
               }
            }
         } else {
            // Use summation over grid points
            const int meshSize = domain_.mesh().size();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < meshSize; ++k) {
                  fExt_ += h_.rgrid(i)[k] * c_.rgrid(i)[k];
               }
            }
            fExt_ /= double(meshSize);
         }
         fExt_ /= phiTot;
         fHelmholtz_ += fExt_;
      }

      // Compute excess interaction free energy [ phi^{T}*chi*phi/2 ]
      if (w_.isSymmetric()) {
         UTIL_CHECK(c_.isAllocatedBasis());
         const int nBasis = domain_.basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (std::abs(chi) > 1.0E-9) {
                  double temp = 0.0;
                  for (int k = 0; k < nBasis; ++k) {
                     temp += c_.basis(i)[k] * c_.basis(j)[k];
                  }
                  if (i == j) {
                     fInter_ += 0.5*chi*temp;
                  } else {
                     fInter_ += chi*temp;
                  }
               }
            }
         }
      } else {
         const int meshSize = domain_.mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (std::abs(chi) > 1.0E-9) {
                  double temp = 0.0;
                  for (int k = 0; k < meshSize; ++k) {
                     temp += c_.rgrid(i)[k] * c_.rgrid(j)[k];
                  }
                  if (i == j) {
                     fInter_ += 0.5*chi*temp;
                  } else {
                     fInter_ += chi*temp;
                  }
               }
            }
         }
         fInter_ /= double(meshSize);
      }
      fInter_ /= phiTot;
      fHelmholtz_ += fInter_;

      // Initialize pressure (-1 x grand-canonical free energy / monomer)
      pressure_ = -fHelmholtz_;

      // Polymer chemical potential corrections to pressure
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            if (PolymerModel::isThread()) {
               length = polymerPtr->length();
            } else {
               length = (double) polymerPtr->nBead();
            }
            if (phi > 1.0E-08) {
               pressure_ += mu * phi / length;
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
            if (phi > 1.0E-08) {
               pressure_ += mu * phi / size;
            }
         }
      }

      hasFreeEnergy_ = true;
   }

   /*
   * Compute SCFT stress for current fields.
   */
   template <int D>
   void System<D>::computeStress()
   {
      if (hasStress_) return;

      stress_.clear();

      // Compute and store stress contribution from Mixture
      if (!mixture_.hasStress()) {
         mixture_.computeStress(mask().phiTot());
      }
      for (int i = 0; i < domain_.unitCell().nParameter(); ++i) {
         stress_.append(mixture_.stress(i));
      }

      if (hasEnvironment()) {
         UTIL_CHECK(!environment().needsUpdate());
         for (int i = 0; i < domain_.unitCell().nParameter(); ++i) {
            if (iterator().flexibleParams()[i]) {
               // Add stress contributions from Environment
               stress_[i] += environment().stress(i);

               // Allow Environment to modify stress contributions to
               // minimize something other than fHelmholtz if desired
               stress_[i] = environment().modifyStress(i,stress_[i]);
            }
         }
      }

      hasStress_ = true;
   }

   /*
   * Write parameter file for SCFT, omitting any sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture_.writeParam(out);
      interaction().writeParam(out);
      domain_.writeParam(out);
      if (hasEnvironment()) {
         environment().writeParam(out);
      }
      if (hasIterator()) {
         iterator().writeParam(out);
      }
      out << "}" << std::endl;
   }

   /*
   * Write thermodynamic properties to file.
   */
   template <int D>
   void System<D>::writeThermo(std::ostream& out)
   {
      if (!hasFreeEnergy_) {
         computeFreeEnergy();
      }

      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
      if (hasExternalFields()) {
         out << "fExt          " << Dbl(fExt_, 18, 11) << std::endl;
      }
      out << std::endl;

      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      if (np > 0) {
         out << "polymers:" << std::endl;
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
         out << "solvents:" << std::endl;
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

      out << "cellParams:" << std::endl;
      for (int i = 0; i < domain_.unitCell().nParameter(); ++i) {
         out << Int(i, 5)
             << "  "
             << Dbl(domain_.unitCell().parameter(i), 18, 11)
             << std::endl;
      }
      out << std::endl;
   }

   /*
   * Write stress properties to file.
   */
   template <int D>
   void System<D>::writeStress(std::ostream& out)
   {
      UTIL_CHECK(mixture_.hasStress()); // ensure stress has been calculated

      out << "stress:";
      if (hasEnvironment()) {
         out << " (Environment contributions not included)";
      }
      out << std::endl;
      
      for (int i = 0; i < domain_.unitCell().nParameter(); ++i) {
         out << Int(i, 5)
            << "  "
            << Dbl(mixture_.stress(i), 18, 11)
            << std::endl;
      }
      out << std::endl;
   }

   // Field Output

   /*
   * Write w fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeWBasis(std::string const & filename) const
   {
      // Preconditions
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());

      domain_.fieldIo().writeFieldsBasis(filename, w_.basis(),
                                         domain_.unitCell());
   }

   /*
   * Write w fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(std::string const & filename) const
   {
      // Preconditions
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(w_.hasData());

      domain_.fieldIo().writeFieldsRGrid(filename, w_.rgrid(),
                                         domain_.unitCell(),
                                         w_.isSymmetric());
   }

   /*
   * Write concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(std::string const & filename) const
   {
      // Preconditions
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedBasis_);
      UTIL_CHECK(c_.hasData());
      UTIL_CHECK(c_.hasData() == hasCFields_);
      UTIL_CHECK(w_.isSymmetric());

      domain_.fieldIo().writeFieldsBasis(filename, c_.basis(),
                                         domain_.unitCell());
   }

   /*
   * Write concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(std::string const & filename) const
   {
      // Preconditions
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(c_.hasData());
      UTIL_CHECK(c_.hasData() == hasCFields_);

      domain_.fieldIo().writeFieldsRGrid(filename, c_.rgrid(),
                                         domain_.unitCell(),
                                         w_.isSymmetric());
   }

   /*
   * Write all concentration fields in real space (r-grid) format, for
   * each block and solvent species, rather than for each monomer type.
   */
   template <int D>
   void System<D>::writeBlockCRGrid(std::string const & filename) const
   {
      // Preconditions
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(c_.hasData());
      UTIL_CHECK(c_.hasData() == hasCFields_);

      // Create array to hold block and solvent c field data
      DArray< RField<D> > blockCFields;
      blockCFields.allocate(mixture_.nSolvent() + mixture_.nBlock());
      int n = blockCFields.capacity();
      for (int i = 0; i < n; i++) {
         blockCFields[i].allocate(domain_.mesh().dimensions());
      }

      // Get c field data from the Mixture
      mixture_.createBlockCRGrid(blockCFields);

      // Write block and solvent c field data to file
      domain_.fieldIo().writeFieldsRGrid(filename, blockCFields,
                                         domain_.unitCell(),
                                         w_.isSymmetric());
   }

   // Propagator Output

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQSlice(std::string const & filename,
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
      UTIL_CHECK(domain_.unitCell().isInitialized());
      Propagator<D> const &
           propagator = polymer.propagator(blockId, directionId);
      RField<D> const & field = propagator.q(segmentId);
      domain_.fieldIo().writeFieldRGrid(filename, field,
                                        domain_.unitCell(),
                                        w_.isSymmetric());
   }

   /*
   * Write the last (tail) slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQTail(std::string const & filename,
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
      UTIL_CHECK(domain_.unitCell().isInitialized());
      RField<D> const &
            field = polymer.propagator(blockId, directionId).tail();
      domain_.fieldIo().writeFieldRGrid(filename, field,
                                        domain_.unitCell(),
                                        w_.isSymmetric());
   }

   /*
   * Write the propagator for a block and direction.
   */
   template <int D>
   void System<D>::writeQ(std::string const & filename,
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
      UTIL_CHECK(domain_.unitCell().isInitialized());
      Propagator<D> const&
           propagator = polymer.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);

      // Write header
      domain_.fieldIo().writeFieldHeader(file, 1, domain_.unitCell(),
                                         w_.isSymmetric());
      file << "mesh" << std::endl
           << "          " << domain_.mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
         file << "slice " << i << std::endl;
         domain_.fieldIo().writeFieldRGrid(file, propagator.q(i),
                                           domain_.unitCell(),
                                           hasHeader);
      }
      file.close();
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
               filename += ".rq";
               writeQ(filename, ip, ib, id);
            }
         }
      }
   }

   // Timer Operations

   /*
   * Write timer values to output stream (computational cost).
   */
   template <int D>
   void System<D>::writeTimers(std::ostream& out)
   {
      if (hasIterator()) {
         iterator().outputTimers(Log::file());
         iterator().outputTimers(out);
      }
      if (hasSimulator()){
         simulator().outputTimers(Log::file());
         simulator().outputTimers(out);
      }
   }

   /*
   * Clear state of all timers.
   */
   template <int D>
   void System<D>::clearTimers()
   {
      if (hasIterator()) {
         iterator().clearTimers();
      }
      if (hasSimulator()){
         simulator().clearTimers();
      }
   }

   // Miscellaneous public functions

   /*
   * Mark c-fields and free energy as outdated/invalid.
   */
   template <int D>
   void System<D>::clearCFields()
   {
      c_.setHasData(false);
      hasCFields_ = false;
      hasFreeEnergy_ = false;
      hasStress_ = false;
   }

   // Private member functions

   /*
   * Allocate memory for fields in grid format.
   */
   template <int D>
   void System<D>::allocateFieldsGrid()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      int nMonomer = mixture_.nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(!isAllocatedGrid_);

      // Alias for mesh dimensions
      IntVec<D> const & dimensions = domain_.mesh().dimensions();

      // Allocate w (chemical potential) fields
      w_.setNMonomer(nMonomer);
      w_.allocateRGrid(dimensions);

      // Allocate c (monomer concentration) fields
      c_.setNMonomer(nMonomer);
      c_.allocateRGrid(dimensions);

      h_.setNMonomer(nMonomer);

      // If hasEnvironment(), allocate mask and h fields
      if (hasEnvironment()) {
         if (environment().generatesMask()) {
            mask_.allocateRGrid(dimensions);
         }
         if (environment().generatesExternalFields()) {
            h_.allocateRGrid(dimensions);
         }
      }

      isAllocatedGrid_ = true;
   }

   /*
   * Allocate memory for fields in basis format.
   */
   template <int D>
   void System<D>::allocateFieldsBasis()
   {
      // Preconditions and constants
      UTIL_CHECK(hasMixture_);
      const int nMonomer = mixture_.nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(domain_.hasGroup());
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = domain_.basis().nBasis();
      UTIL_CHECK(nBasis > 0);
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(!isAllocatedBasis_);

      // Allocate basis fields in w and c field containers
      w_.allocateBasis(nBasis);
      c_.allocateBasis(nBasis);

      // If hasEnvironment(), allocate mask and h fields
      if (hasEnvironment()) {
         if (environment().generatesMask()) {
            mask_.allocateBasis(nBasis);
         }
         if (environment().generatesExternalFields()) {
            h_.allocateBasis(nBasis);
         }
      }

      isAllocatedBasis_ = true;
   }

   /*
   * Read a filename string and echo to log file (used in readCommands).
   */
   template <int D>
   void System<D>::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      if (in.fail()) {
          UTIL_THROW("Unable to read a string parameter.");
      }
      Log::file() << " " << Str(string, 20) << std::endl;
   }

   /*
   * Read floating point number, echo to log file (used in readCommands).
   */
   template <int D>
   void System<D>::readEcho(std::istream& in, double& value) const
   {
      in >> value;
      if (in.fail()) {
          UTIL_THROW("Unable to read floating point parameter.");
      }
      Log::file() << " " << Dbl(value, 20) << std::endl;
   }

} // namespace Rpc
} // namespace Pscf
#endif
