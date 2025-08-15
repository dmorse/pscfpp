#ifndef RPC_SYSTEM_TPP
#define RPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field
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
#include <rpc/scft/ScftThermo.h>
#include <rpc/environment/EnvironmentFactory.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>


#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>
#include <prdc/environment/Environment.h>

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
      mixtureModifierPtr_(nullptr),
      interactionPtr_(nullptr),
      environmentPtr_(nullptr),
      environmentFactoryPtr_(nullptr),
      scftPtr_(nullptr),
      iteratorPtr_(nullptr),
      iteratorFactoryPtr_(nullptr),
      sweepPtr_(nullptr),
      sweepFactoryPtr_(nullptr),
      simulatorPtr_(nullptr),
      simulatorFactoryPtr_(nullptr),
      polymerModel_(PolymerModel::Thread),
      isAllocatedGrid_(false),
      isAllocatedBasis_(false),
      hasMixture_(false)
   {
      setClassName("System");  // Set block label in parameter file
      BracketPolicy::set(BracketPolicy::Optional);


      // Create dynamically allocated objects owned by this System
      mixtureModifierPtr_ = new MixtureModifier<D>();
      interactionPtr_ = new Interaction();
      environmentFactoryPtr_ = new EnvironmentFactory<D>(*this);
      scftPtr_ = new ScftThermo<D>(*this);
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this);
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      simulatorFactoryPtr_ = new SimulatorFactory<D>(*this);

      // Create associations among class members
      mixtureModifier().associate(mixture_);
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

      // Note the arguments of setReadUnitCell functions used above:
      // When w_ is read from a file in basis or r-grid format, the
      // parameters of the system unit cell, domain_.unitCell(), are
      // reset to those in the field file header. When imposed fields h_
      // and mask_ are read from a file, however, unit cell parameters
      // from the field file header are read into a mutable workspace
      // object, tmpUnitCell_, and ultimately discarded.

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
      if (scftPtr_) {
         delete scftPtr_;
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
      mixture_.setFieldIo(domain_.fieldIo());
      mixture_.allocate();

      // Allocate memory for field containers in r-grid form
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
      if (environmentPtr_) {
         environment().setNParams(domain_.unitCell().nParameter());
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
            w().readBasis(filename);
            UTIL_CHECK(domain_.unitCell().isInitialized());
            UTIL_CHECK(domain_.basis().isInitialized());
            UTIL_CHECK(!domain_.waveList().hasKSq());
            UTIL_CHECK(isAllocatedBasis_);
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(!scft().hasData());
         } else
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);
            w().readRGrid(filename);
            UTIL_CHECK(domain_.unitCell().isInitialized());
            UTIL_CHECK(!domain_.waveList().hasKSq());
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(!scft().hasData());
         } else
         if (command == "SET_UNIT_CELL") {
            UnitCell<D> unitCell;
            in >> unitCell;
            Log::file() << "   " << unitCell << std::endl;
            setUnitCell(unitCell);
            UTIL_CHECK(domain_.unitCell().isInitialized());
            UTIL_CHECK(!domain_.waveList().hasKSq());
            UTIL_CHECK(!c_.hasData());
            UTIL_CHECK(!scft().hasData());
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
            scft().write(file);
            file.close();
         } else
         if (command == "WRITE_STRESS") {
            readEcho(in, filename);
            std::ofstream file;
            if (!mixture().hasStress()) {
               computeStress();
            }
            fileMaster_.openOutputFile(filename, file,
                                       std::ios_base::app);
            mixture().writeStress(file);
            if (hasEnvironment()) {
               environment().writeStress(file);
            }
            file.close();
         } else
         if (command == "WRITE_W_BASIS") {
            readEcho(in, filename);
            w().writeBasis(filename);
         } else
         if (command == "WRITE_W_RGRID") {
            readEcho(in, filename);
            w().writeRGrid(filename);
         } else
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            c().writeBasis(filename);
         } else
         if (command == "WRITE_C_RGRID") {
            readEcho(in, filename);
            c().writeRGrid(filename);
         } else
         if (command == "WRITE_BLOCK_C_RGRID") {
            readEcho(in, filename);
            mixture_.writeBlockCRGrid(filename);
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
            mixture_.writeQSlice(filename, polymerId, blockId,
                                 directionId, segmentId);
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
            mixture_.writeQTail(filename, polymerId, blockId, directionId);
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
            mixture_.writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(in, filename);
            mixture_.writeQAll(filename);
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
         if (command == "ESTIMATE_W_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo.estimateWBasis(inFileName, outFileName,
                                   interaction().chi());
         } else
         if (command == "READ_H_BASIS") {
            readEcho(in, filename);
            h().readBasis(filename);
            UTIL_CHECK(!c_.hasData());
         } else
         if (command == "READ_H_RGRID") {
            readEcho(in, filename);
            h().readRGrid(filename);
            UTIL_CHECK(!c_.hasData());
         } else
         if (command == "WRITE_H_BASIS") {
            readEcho(in, filename);
            h().writeBasis(filename);
         } else
         if (command == "WRITE_H_RGRID") {
            readEcho(in, filename);
            h().writeRGrid(filename);
         } else
         if (command == "READ_MASK_BASIS") {
            readEcho(in, filename);
            mask().readBasis(filename);
         } else
         if (command == "READ_MASK_RGRID") {
            readEcho(in, filename);
            mask().readRGrid(filename);
         } else
         if (command == "WRITE_MASK_BASIS") {
            readEcho(in, filename);
            mask().writeBasis(filename);
         } else
         if (command == "WRITE_MASK_RGRID") {
            readEcho(in, filename);
            mask().writeRGrid(filename);
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
      mixture_.setIsSymmetric(w_.isSymmetric());
      c_.setHasData(true);
      scft().clear();

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
   * Compute SCFT stress for current fields.
   */
   template <int D>
   void System<D>::computeStress()
   {
      // Compute and store standard Mixture stress
      if (!mixture_.hasStress()) {
         mixture_.computeStress(mask().phiTot());
      }

      // If necessary, compute and store Environment stress
      if (hasEnvironment()) {
         if (!environment().hasStress()) {
            environment().computeStress(mixture(),
                                        iterator().flexibleParams());
         }
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

      // If converged, compute related thermodynamic properties
      if (!error) {
         scft().compute();
         scft().write(Log::file());
         if (!iterator().isFlexible()) {
            if (!mixture().hasStress()) {
               computeStress();
            }
            mixture().writeStress(Log::file());
            if (hasEnvironment()) {
               environment().writeStress(Log::file());
            }
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

   /*
   * Mark c-fields and free energy as outdated/invalid.
   */
   template <int D>
   void System<D>::clearCFields()
   {
      c_.setHasData(false);
      scft().clear();
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

   // Property Output

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

   // Timer Operations

   /*
   * Write timer values to output stream (computational cost).
   */
   template <int D>
   void System<D>::writeTimers(std::ostream& out) const
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

      // Allocate c (monomer concentration) fields
      c_.setNMonomer(nMonomer);
      c_.allocateRGrid(dimensions);

      // Allocate w (chemical potential) fields
      w_.setNMonomer(nMonomer);
      w_.allocateRGrid(dimensions);

      // Set nMonomer for external field container
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
