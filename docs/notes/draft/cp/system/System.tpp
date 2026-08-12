#ifndef PRDC_CL_SYSTEM_TPP
#define PRDC_CL_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/math/IntVec.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/param/BracketPolicy.h>
#include <util/param/ParamComponent.h>
#include <util/signal/Signal.h>
#include <util/format/Str.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <util/misc/FileMaster.h>

#include <string>
#include <unistd.h>

namespace Pscf {
namespace Cp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   System<D,T>::System(typename T::System& system)
    : w_(),
      c_(),
      systemPtr_(&system),
      mixturePtr_(nullptr),
      //mixtureModifierPtr_(nullptr),
      domainPtr_(nullptr),
      interactionPtr_(nullptr),
      //simulatorPtr_(nullptr),
      //simulatorFactoryPtr_(nullptr),
      fileMasterPtr_(nullptr),
      polymerModel_(PolymerModel::Thread),
      isAllocated_(false),
      hasMixture_(false)
   {
      ParamComposite::setClassName("System");  // label in parameter file
      BracketPolicy::set(BracketPolicy::Optional);

      // Create dynamically allocated objects owned by this System
      mixturePtr_ = new typename T::Mixture();
      // mixtureModifierPtr_ = new typename T::MixtureModifier();
      interactionPtr_ = new typename T::Interaction();
      domainPtr_ = new typename T::Domain();
      //simulatorFactoryPtr_ = new typename T::SimulatorFactory(*systemPtr_);
      fileMasterPtr_ = new FileMaster();

      // Create associations among class members
      //mixtureModifier().associate(mixture_());
      domain_().setFileMaster(*fileMasterPtr_);
      w_.setFieldIo(domain_().fieldIo());
      w_.setReadUnitCell(domain_().unitCell());
      w_.setWriteUnitCell(domain_().unitCell());
      c_.setFieldIo(domain_().fieldIo());
      c_.setWriteUnitCell(domain_().unitCell());

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
      Signal<void>& cellSignal = domain_().unitCell().signal();
      cellSignal.addObserver(*this, &System<D,T>::clearUnitCellData);

      // Signal triggered by w-field modification
      w_.signal().addObserver(*this, &System<D,T>::clearCFields);

   }

   /*
   * Destructor.
   */
   template <int D, class T>
   System<D,T>::~System()
   {
      delete mixturePtr_;
      //delete mixtureModifierPtr_;
      delete interactionPtr_;
      delete domainPtr_;
      //delete simulatorFactoryPtr_;
      //if (simulatorPtr_) {
      //   delete simulatorPtr_;
      //}
      delete fileMasterPtr_;
   }

   // Lifetime (called in main program)

   /*
   * Process command line options.
   */
   template <int D, class T>
   void System<D,T>::setOptions(int argc, char **argv)
   {
      bool eFlag = false;  // echo
      bool dFlag = false;  // spatial dimension (1, 2, or 3)
      bool pFlag = false;  // param file
      bool cFlag = false;  // command file
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      bool tFlag = false;  // thread count
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
         case 't': // thread count, if set by user
            tFlag = true;
            tArg = atoi(optarg);
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }

      // Process program arguments

      // Check required option -d
      if (dFlag) {
         UTIL_CHECK(D == dArg);
      } else {
         UTIL_THROW("Missing required -d option");
      }

      // Check required option -p, set parameter file name
      if (pFlag) {
         fileMaster().setParamFileName(std::string(pArg));
      } else {
         UTIL_THROW("Missing required -p option - no parameter file");
      }

      // Check required option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      } else {
         UTIL_THROW("Missing required -c option - no command file");
      }

      // If option -e, set to echo parameters as they are read
      if (eFlag) {
         Util::ParamComponent::setEcho(true);
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

      // If option -t, process the thread count
      if (tFlag) {
         if (tArg <= 0) {
            UTIL_THROW("Error: Non-positive thread count -t option");
         }
         setThreadCount(tArg);  // Default implementation does nothing
      }
   }

   /*
   * Read parameters and initialize.
   */
   template <int D, class T>
   void System<D,T>::readParameters(std::istream& in)
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
      // that it was not reset within the remainder of this function.

      int nSetPolymerModel = PolymerModel::nSet();

      // Read the Mixture{ ... } block
      readParamComposite(in, mixture_());
      hasMixture_ = true;

      int nm = mixture_().nMonomer();   // number of monomer types
      int np = mixture_().nPolymer();   // number of polymer species
      int ns = mixture_().nSolvent();   // number of solvent species
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(np > 0);
      UTIL_CHECK(ns >= 0);
      UTIL_CHECK(np + ns > 0);

      // Read the Interaction{ ... } block
      interaction().setNMonomer(nm);
      readParamComposite(in, interaction());

      // Read the Domain{ ... } block
      readParamComposite(in, domain_());
      UTIL_CHECK(domain_().mesh().size() > 0);
      UTIL_CHECK(domain_().unitCell().nParameter() > 0);
      UTIL_CHECK(domain_().unitCell().lattice() != UnitCell<D>::Null);
      domain_().fieldIo().setNMonomer(nm);

      // Complete setup and allocation of the Mixture
      mixture_().associate(domain_().mesh(), domain_().fft(),
                           domain_().unitCell(), domain_().waveList());
      mixture_().setFieldIo(domain_().fieldIo());
      mixture_().allocate();

      // Allocate memory for w- and c-field containers 
      allocateFields();

      #if 0
      // Optionally read and construct a Simulator 
      if (!isEnd) {
         simulatorPtr_ =
            simulatorFactoryPtr_->readObjectOptional(in, *this,
                                                     className, isEnd);
         if (!simulatorPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Simulator{ [absent] }\n";
         }
      }
      #endif

      // Check that the polymer model was not modified after initialization
      UTIL_CHECK(PolymerModel::nSet() == nSetPolymerModel);

      // Note: The main pscf_cpc or pscf_cpg program can irreversibly lock the
      // global polymer model after reading the parameter file, to prohibit 
      // later changes during processing of the command file. The polymer model
      // enumeration is not locked here so that this function can be called 
      // repeatedly during unit testing.
   }

   /*
   * Read parameter file (including open and closing brackets).
   */
   template <int D, class T>
   void System<D,T>::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());
      readParameters(in);
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   template <int D, class T>
   void System<D,T>::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read and execute commands from a specified command file.
   */
   template <int D, class T>
   void System<D,T>::readCommands(std::istream &in)
   {
      UTIL_CHECK(isAllocated_);
      std::string command, filename, inFileName, outFileName;
      //typename T::FieldIo const & fieldIo = domain_().fieldIo();

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
         if (command == "READ_W") {
            readEcho(in, filename);
            w().readFields(filename);
            UTIL_CHECK(domain_().unitCell().isInitialized());
            UTIL_CHECK(!domain_().waveList().hasKSq());
            UTIL_CHECK(!c_.hasData());
         } else
         if (command == "SET_UNIT_CELL") {
            UnitCell<D> unitCell;
            in >> unitCell;
            Log::file() << "   " << unitCell << std::endl;
            setUnitCell(unitCell);
         } else
         if (command == "COMPUTE") {
            // Solve the modified diffusion equation for the mixture
            compute();
         } else

         #if 0
         if (command == "SIMULATE") {
            // Perform a field theoretic simulation of nStep steps
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
         #endif

         if (command == "WRITE_PARAM") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file);
            writeParam(file);
            file.close();
         } else
         if (command == "WRITE_W") {
            readEcho(in, filename);
            w().writeFields(filename);
         } else
         if (command == "WRITE_C") {
            readEcho(in, filename);
            c().writeFields(filename);
         } else

         #if 0
         if (command == "WRITE_BLOCK_C_RGRID") {
            readEcho(in, filename);
            mixture_().writeBlockCRGrid(filename);
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
            mixture_().writeQSlice(filename, polymerId, blockId,
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
            mixture_().writeQTail(filename, polymerId, blockId, directionId);
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
            mixture_().writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(in, filename);
            mixture_().writeQAll(filename);
         } else
         #endif

         #if 0
         if (command == "WRITE_TIMERS") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file);
            writeTimers(file);
            file.close();
         } else
         if (command == "CLEAR_TIMERS") {
            clearTimers();
         } else 
         #endif

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
   template <int D, class T>
   void System<D,T>::readCommands()
   {
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile());
   }

   // Primary Field Theory Computations

   /*
   * Solve MDE for current w fields, without iteration.
   */
   template <int D, class T>
   void System<D,T>::compute(bool needStress)
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(w_.isAllocated());
      UTIL_CHECK(c_.isAllocated());
      UTIL_CHECK(w_.hasData());
      clearCFields();

      // Solve the modified diffusion equation (without iteration)
      mixture_().compute(w_.fields(), c_.fields());
      c_.setHasData(true);

   }

   #if 0
   /*
   * Perform a field theoretic simulation of nStep steps.
   */
   template <int D, class T>
   void System<D,T>::simulate(int nStep)
   {
      UTIL_CHECK(nStep > 0);
      UTIL_CHECK(hasSimulator());
      clearCFields();

      simulator().simulate(nStep);
   }
   #endif

   /*
   * Mark c-fields and free energy as outdated/invalid.
   */
   template <int D, class T>
   void System<D,T>::clearCFields()
   {  c_.setHasData(false); }

   // Unit Cell Modifiers

   /*
   * Set the system unit cell.
   */
   template <int D, class T>
   void System<D,T>::setUnitCell(UnitCell<D> const & unitCell)
   {
      // Preconditions
      UTIL_CHECK(domain_().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(domain_().lattice() == unitCell.lattice());

      // Set system unit cell, using UnitCell<D> assignment operator
      domain_().unitCell() = unitCell;

      // Postconditions
      UTIL_CHECK(domain_().unitCell().isInitialized());
      UTIL_CHECK(domain_().unitCell().lattice() == domain_().lattice());
      UTIL_CHECK(!domain_().waveList().hasKSq());
   }

   /*
   * Set parameters of the system unit cell.
   */
   template <int D, class T>
   void System<D,T>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      // Precondition
      UTIL_CHECK(domain_().lattice() != UnitCell<D>::Null);

      // Set system unit cell
      if (domain_().unitCell().lattice() == UnitCell<D>::Null) {
         domain_().unitCell().set(domain_().lattice(), parameters);
      } else {
         domain_().unitCell().setParameters(parameters);
      }

      // Postconditions
      UTIL_CHECK(domain_().unitCell().isInitialized());
      UTIL_CHECK(domain_().unitCell().lattice() == domain_().lattice());
      UTIL_CHECK(!domain_().waveList().hasKSq());
   }

   /*
   * Notify System members that unit cell parameters have been modified.
   */
   template <int D, class T>
   void System<D,T>::clearUnitCellData()
   {
      clearCFields();
      mixture_().clearUnitCellData();
      domain_().waveList().clearUnitCellData();
   }

   // Timer Operations

   #if 0
   /*
   * Write timer values to output stream (computational cost).
   */
   template <int D, class T>
   void System<D,T>::writeTimers(std::ostream& out) const
   {
      if (hasSimulator()){
         simulator().outputTimers(Log::file());
         simulator().outputTimers(out);
      }
   }

   /*
   * Clear state of all timers.
   */
   template <int D, class T>
   void System<D,T>::clearTimers()
   {
      if (hasSimulator()){
         simulator().clearTimers();
      }
   }
   #endif

   // Private member functions

   /*
   * Allocate memory for fields in grid format.
   */
   template <int D, class T>
   void System<D,T>::allocateFields()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      int nMonomer = mixture_().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(domain_().mesh().size() > 0);
      UTIL_CHECK(!isAllocated_);

      // Alias for mesh dimensions
      IntVec<D> const & dimensions = domain_().mesh().dimensions();

      // Allocate w (chemical potential) fields
      //w_.setNMonomer(nMonomer);
      w_.allocate(nMonomer, dimensions);

      // Allocate c (monomer concentration) fields
      //c_.setNMonomer(nMonomer);
      c_.allocate(nMonomer, dimensions);

      isAllocated_ = true;
   }

   /*
   * Read a filename string and echo to log file (used in readCommands).
   */
   template <int D, class T>
   void System<D,T>::readEcho(std::istream& in, std::string& string) const
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
   template <int D, class T>
   void System<D,T>::readEcho(std::istream& in, double& value) const
   {
      in >> value;
      if (in.fail()) {
          UTIL_THROW("Unable to read floating point parameter.");
      }
      Log::file() << " " << Dbl(value, 20) << std::endl;
   }

} // namespace Cp
} // namespace Pscf
#endif
