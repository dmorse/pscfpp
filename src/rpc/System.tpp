#ifndef RPC_SYSTEM_TPP
#define RPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <rpc/simulate/Simulator.h>
#include <rpc/simulate/SimulatorFactory.h>
#include <rpc/simulate/compressor/Compressor.h>
#include <rpc/sweep/Sweep.h>
#include <rpc/sweep/SweepFactory.h>
#include <rpc/iterator/Iterator.h>
#include <rpc/iterator/IteratorFactory.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>

#include <pscf/inter/Interaction.h>
#include <pscf/math/IntVec.h>
#include <pscf/homogeneous/Clump.h>

#include <util/param/BracketPolicy.h>
#include <util/param/ParamComponent.h>
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
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      simulatorPtr_(0),
      simulatorFactoryPtr_(0),
      w_(),
      c_(),
      h_(),
      mask_(),
      fHelmholtz_(0.0),
      fIdeal_(0.0),
      fInter_(0.0),
      fExt_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocatedGrid_(false),
      isAllocatedBasis_(false),
      hasCFields_(false),
      hasFreeEnergy_(false)
   {
      setClassName("System");
      domain_.setFileMaster(fileMaster_);
      w_.setFieldIo(domain_.fieldIo());
      h_.setFieldIo(domain_.fieldIo());
      mask_.setFieldIo(domain_.fieldIo());
      interactionPtr_ = new Interaction();
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this);
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      simulatorFactoryPtr_ = new SimulatorFactory<D>(*this);
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
      if (simulatorPtr_) {
         delete simulatorPtr_;
      }
      if (simulatorFactoryPtr_) {
         delete simulatorFactoryPtr_;
      }
   }

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
         case 't':
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

      // Check existence and consistency of -d option
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
         } else {
            // Set thread count
         }
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

      // Read the Interaction{ ... } block
      interaction().setNMonomer(nm);
      readParamComposite(in, interaction());

      // Read the Domain{ ... } block
      readParamComposite(in, domain_);
      UTIL_CHECK(domain().mesh().size() > 0);
      UTIL_CHECK(domain().unitCell().nParameter() > 0);
      UTIL_CHECK(domain().unitCell().lattice() != UnitCell<D>::Null);

      // Setup mixture
      mixture_.setDiscretization(domain().mesh(), domain().fft());
      mixture_.setupUnitCell(domain().unitCell());

      // Allocate memory for w and c fields
      allocateFieldsGrid();
      if (domain().basis().isInitialized()) {
         allocateFieldsBasis();
      }

      // Optionally instantiate an Iterator object
      std::string className;
      bool isEnd;
      iteratorPtr_ =
         iteratorFactoryPtr_->readObjectOptional(in, *this, className,
                                                 isEnd);
      if (!iteratorPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Iterator{ [absent] }\n";
      }

      // Optionally instantiate a Sweep object (if an iterator exists)
      if (iteratorPtr_ && !isEnd) {
         sweepPtr_ =
            sweepFactoryPtr_->readObjectOptional(in, *this, className,
                                                 isEnd);
         if (!sweepPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Sweep{ [absent] }\n";
         }
      }

      // Optionally instantiate a Simulator object
      if (!isEnd) {
         simulatorPtr_ =
            simulatorFactoryPtr_->readObjectOptional(in, *this,
                                                     className, isEnd);
         if (!simulatorPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Simulator{ [absent] }\n";
         }
      }

      // Initialize homogeneous object
      // NOTE: THIS OBJECT IS NOT USED AT ALL.
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

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
            readWBasis(filename);
         } else
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);
            readWRGrid(filename);
         } else
         if (command == "ESTIMATE_W_FROM_C") {
            readEcho(in, inFileName);
            estimateWfromC(inFileName);
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
         if (command == "WRITE_TIMERS") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster_.openOutputFile(filename, file);
            writeTimers(file);
            file.close();
         } else
         if (command == "CLEAR_TIMERS") {
            clearTimers();
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
                        << Str("direction ID  ", 21) << directionId << "\n"
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
            writeStars(filename);
         } else
         if (command == "WRITE_WAVES") {
            readEcho(in, filename);
            writeWaves(filename);
         } else
         if (command == "WRITE_GROUP") {
            readEcho(in, filename);
            writeGroup(filename);
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
            double epsilon;
            readEcho(in, inFileName);
            readEcho(in, epsilon);
            bool hasSymmetry;
            hasSymmetry = checkRGridFieldSymmetry(inFileName, epsilon);
            if (hasSymmetry) {
               Log::file() << std::endl
                   << "Symmetry of r-grid file matches this space group."
                   << std::endl << std::endl;
            } else {
               Log::file() << std::endl
                 << "Symmetry of r-grid file does not match this\n"
                 << "space group to within error threshold of "
                 << Dbl(epsilon) << " ." << std::endl << std::endl;
            }
         } else
         if (command == "COMPARE_BASIS") {

            // Get two filenames for comparison
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);

            DArray< DArray<double> > Bfield1, Bfield2;
            UnitCell<D> tmpUnitCell;
            domain().fieldIo().readFieldsBasis(filecompare1, Bfield1,
                                               tmpUnitCell);
            domain().fieldIo().readFieldsBasis(filecompare2, Bfield2,
                                               tmpUnitCell);
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
            UnitCell<D> tmpUnitCell;
            domain_.fieldIo().readFieldsRGrid(filecompare1, Rfield1,
                                              tmpUnitCell);
            domain_.fieldIo().readFieldsRGrid(filecompare2, Rfield2,
                                              tmpUnitCell);
            // Note: Rfield1, Rfield2 will be allocated by readFieldsRGrid

            // Compare and output report
            compare(Rfield1, Rfield2);

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

            expandRGridDimension(inFileName, outFileName,
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

            replicateUnitCell(inFileName, outFileName, replicas);
         } else
         if (command == "READ_H_BASIS") {
            readEcho(in, filename);
            if (!h_.isAllocatedBasis()) {
               h_.allocateBasis(domain().basis().nBasis());
            }
            if (!h_.isAllocatedRGrid()) {
               h_.allocateRGrid(domain().mesh().dimensions());
            }
            UnitCell<D> tmpUnitCell;
            h_.readBasis(filename, tmpUnitCell);
         } else
         if (command == "READ_H_RGRID") {
            readEcho(in, filename);
            if (!h_.isAllocatedRGrid()) {
               h_.allocateRGrid(domain().mesh().dimensions());
            }
            UnitCell<D> tmpUnitCell;
            h_.readRGrid(filename, tmpUnitCell);
         } else
         if (command == "WRITE_H_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            UTIL_CHECK(h_.isSymmetric());
            domain().fieldIo().writeFieldsBasis(filename, h_.basis(), 
                                                domain().unitCell());
         } else
         if (command == "WRITE_H_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            domain().fieldIo().writeFieldsRGrid(filename, h_.rgrid(), 
                                                domain().unitCell());
         } else
         if (command == "READ_MASK_BASIS") {
            UTIL_CHECK(domain_.basis().isInitialized());
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(domain().basis().nBasis(), 
                              domain().mesh().dimensions());
            }
            UnitCell<D> tmpUnitCell;
            mask_.readBasis(filename, tmpUnitCell);
         } else
         if (command == "READ_MASK_RGRID") {
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(domain().basis().nBasis(), 
                              domain().mesh().dimensions());
            }
            UnitCell<D> tmpUnitCell;
            mask_.readBasis(filename, tmpUnitCell);
         } else
         if (command == "WRITE_MASK_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            UTIL_CHECK(mask_.isSymmetric());
            domain().fieldIo().writeFieldBasis(filename, mask_.basis(), 
                                               domain().unitCell());
         } else
         if (command == "WRITE_MASK_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            domain().fieldIo().writeFieldRGrid(filename, mask_.rgrid(), 
                                               domain().unitCell(),
                                      mask_.isSymmetric());
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

   // W Field Modifier Functions

   /*
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      // Precondition
      UTIL_CHECK(domain_.hasGroup());

      if (!domain().unitCell().isInitialized()) {
         readFieldHeader(filename);
      }
      UTIL_CHECK(domain().unitCell().isInitialized());
      UTIL_CHECK(domain().basis().isInitialized());
      UTIL_CHECK(domain().basis().nBasis() > 0);
      if (!isAllocatedBasis_) {
         allocateFieldsBasis();
      }

      // Read w fields
      w_.readBasis(filename, domain_.unitCell());
      mixture_.setupUnitCell(domain().unitCell());

      hasCFields_ = false;
      hasFreeEnergy_ = false;

   }

   /*
   * Read w-fields in real-space grid (r-grid) format.
   */
   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {
      UTIL_CHECK(isAllocatedGrid_);

      // If necessary, peek at header to initialize unit cell
      if (!domain().unitCell().isInitialized()) {
         readFieldHeader(filename);
      }
      UTIL_CHECK(domain().unitCell().isInitialized());
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         UTIL_CHECK(domain_.basis().isInitialized());
         UTIL_CHECK(domain_.basis().nBasis() > 0);
         allocateFieldsBasis();
      }

      // Read w fields
      w_.readRGrid(filename, domain_.unitCell());

      // Update UnitCell in Mixture
      mixture_.setupUnitCell(domain().unitCell());

      hasCFields_ = false;
      hasFreeEnergy_ = false;
   }

   /*
   * Construct estimate for w fields from c fields, by setting xi=0.
   *
   * Modifies wFields and wFieldsRGrid.
   */
   template <int D>
   void System<D>::estimateWfromC(std::string const & filename)
   {
      UTIL_CHECK(hasMixture_);
      const int nm = mixture_.nMonomer();
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(domain_.hasGroup());

      if (!domain().unitCell().isInitialized()) {
         readFieldHeader(filename);
      }
      UTIL_CHECK(domain().unitCell().isInitialized());
      UTIL_CHECK(domain_.basis().isInitialized());
      if (!isAllocatedBasis_) {
         allocateFieldsBasis();
      }
      const int nb = domain_.basis().nBasis();
      UTIL_CHECK(nb > 0);

      // Read c fields into temporary array and set unit cell
      domain_.fieldIo().readFieldsBasis(filename, tmpFieldsBasis_,
                                        domain_.unitCell());

      // Update UnitCell in Mixture
      mixture_.setupUnitCell(domain().unitCell());

      DArray<double> wtmp;
      wtmp.allocate(nm);

      // Compute estimated w fields from c fields
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

      // Store initial guess for w fields
      w_.setBasis(tmpFieldsBasis_);

      hasCFields_ = false;
      hasFreeEnergy_ = false;

   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      UTIL_CHECK(domain_.hasGroup());
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(isAllocatedBasis_);
      w_.setBasis(fields);
      hasCFields_ = false;
      hasFreeEnergy_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray< RField<D> > const & fields)
   {
      UTIL_CHECK(isAllocatedGrid_);
      w_.setRGrid(fields);
      hasCFields_ = false;
      hasFreeEnergy_ = false;
   }

   // Unit Cell Modifiers

   /*
   * Set the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      domain_.setUnitCell(unitCell);
      mixture_.setupUnitCell(domain().unitCell());
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         allocateFieldsBasis();
      }
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
      mixture_.setupUnitCell(domain().unitCell());
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         allocateFieldsBasis();
      }
   }

   /*
   * Set parameters of the system unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      domain_.setUnitCell(parameters);
      mixture_.setupUnitCell(domain().unitCell());
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         allocateFieldsBasis();
      }
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
      hasFreeEnergy_ = false;

      // If w fields are symmetric, compute basis components for c-fields
      if (w_.isSymmetric()) {
         UTIL_CHECK(c_.isAllocatedBasis());
         domain_.fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis(), false);
      }

      // Compute stress if needed
      if (needStress) {
         mixture_.computeStress();
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
      if (iterator().isSymmetric()) {
         UTIL_CHECK(w_.isSymmetric());
      }
      hasCFields_ = false;
      hasFreeEnergy_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator (return 0 for convergence, 1 for failure)
      int error = iterator().solve(isContinuation);
      hasCFields_ = true;

      // If converged, compute related properties
      if (!error) {
         computeFreeEnergy();  // Sets hasFreeEnergy_ = true
         if (!iterator().isFlexible()) {
            mixture_.computeStress();
         }
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

   /*
   * Perform a field theoretic simulation of nStep steps.
   */
   template <int D>
   void System<D>::simulate(int nStep)
   {
      UTIL_CHECK(nStep > 0);
      UTIL_CHECK(hasSimulator());
      hasCFields_ = false;
      hasFreeEnergy_ = false;

      simulator().simulate(nStep);
      hasCFields_ = true;
   }

   /*
   * Expand the dimensions of RField
   */
   template <int D>
   void System<D>::expandRGridDimension(const std::string & inFileName,
                                        const std::string & outFileName,
                                        int d,
                                        DArray<int> newGridDimensions)
   {
      UTIL_CHECK(d > D);

      // Read fields
      UnitCell<D> tmpUnitCell;
      domain().fieldIo().readFieldsRGrid(inFileName, 
                                         tmpFieldsRGrid_, 
                                         tmpUnitCell);

      // Expand Fields
      domain_.fieldIo().expandRGridDimension(outFileName,
                                              tmpFieldsRGrid_,
                                              tmpUnitCell,
                                              d, newGridDimensions);
   }

   /*
   * Replicate unit cell
   */
   template <int D>
   void System<D>::replicateUnitCell(const std::string & inFileName,
                                     const std::string & outFileName,
                                     IntVec<D> const & replicas)
   {
      // Read fields
      UnitCell<D> tmpUnitCell;
      domain().fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, 
                                         tmpUnitCell);

      // Replicate Fields
      domain_.fieldIo().replicateUnitCell(outFileName,
                                          tmpFieldsRGrid_,
                                          tmpUnitCell,
                                          replicas);
   }

   // Thermodynamic Properties

   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(!hasFreeEnergy_);

      int np = mixture_.nPolymer();   // Number of polymer species
      int ns = mixture_.nSolvent();   // Number of solvent species
      int nm  = mixture_.nMonomer();  // Number of monomer types

      // Initialize all free energy contributions to zero
      fHelmholtz_ = 0.0;
      fIdeal_ = 0.0;
      fInter_ = 0.0;
      fExt_ = 0.0;

      double phi, mu;

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

      // Volume integrals with a mask: If the system has a mask, then the
      // volume that should be used in calculating free energy/pressure
      // is the volume available to the polymers, not the total unit cell
      // volume. We thus divide all terms that involve integrating over
      // the unit cell volume by quantity mask().phiTot(), which is the
      // volume fraction of the unit cell that is occupied by material.
      // This properly scales them to the correct value. fExt_, fInter_,
      // and the Legendre transform component of fIdeal_ all require
      // this scaling. If no mask is present, mask.phiTot() = 1 and no
      // scaling occurs.

      // Compute Legendre transform subtraction from fIdeal_
      double temp = 0.0;
      if (w_.isSymmetric()) {
         // Use expansion in symmetry-adapted orthonormal basis
         const int nBasis = domain_.basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               temp -= w_.basis(i)[k] * c_.basis(i)[k];
            }
         }
      } else {
         // Use summation over grid points
         const int meshSize = domain().mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < meshSize; ++k) {
               temp -= w_.rgrid(i)[k] * c_.rgrid(i)[k];
            }
         }
         temp /= double(meshSize);
      }
      temp /= mask().phiTot();
      fIdeal_ += temp;
      fHelmholtz_ += fIdeal_;

      // Compute contribution from external fields, if they exist
      if (hasExternalFields()) {
         if (w_.isSymmetric()) {
            // Use expansion in symmetry-adapted orthonormal basis
            const int nBasis = domain_.basis().nBasis();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < nBasis; ++k) {
                  fExt_ += h_.basis(i)[k] * c_.basis(i)[k];
               }
            }
         } else {
            // Use summation over grid points
            const int meshSize = domain().mesh().size();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < meshSize; ++k) {
                  fExt_ += h_.rgrid(i)[k] * c_.rgrid(i)[k];
               }
            }
            fExt_ /= double(meshSize);
         }
         fExt_ /= mask().phiTot();
         fHelmholtz_ += fExt_;
      }

      // Compute excess interaction free energy [ phi^{T}*chi*phi/2 ]
      if (w_.isSymmetric()) {
         const int nBasis = domain_.basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (abs(chi) > 1.0E-9) {
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
         const int meshSize = domain().mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (abs(chi) > 1.0E-9) {
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
      fInter_ /= mask().phiTot();
      fHelmholtz_ += fInter_;

      // Initialize pressure
      pressure_ = -fHelmholtz_;

      // Polymer chemical potential corrections to pressure
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            if (phi > 1.0E-08) {
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
            if (phi > 1.0E-08) {
               pressure_ += mu * phi /size;
            }
         }
      }

      hasFreeEnergy_ = true;
   }

   // Output Operations

   /*
   * Write time cost to file.
   */
   template <int D>
   void System<D>::writeTimers(std::ostream& out)
   {
      if (iteratorPtr_) {
         iterator().outputTimers(Log::file());
         iterator().outputTimers(out);
      }
      if (hasSimulator()){
         simulator().outputTimers(Log::file());
         simulator().outputTimers(out);
      }
   }

   /*
   * Clear timers.
   */
   template <int D>
   void System<D>::clearTimers()
   {
      if (iteratorPtr_) {
         iterator().clearTimers();
      }
      if (hasSimulator()){
         simulator().clearTimers();
      }
   }

   /*
   * Write parameter file, omitting any sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture_.writeParam(out);
      interaction().writeParam(out);
      domain_.writeParam(out);
      if (iteratorPtr_) {
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
      for (int i = 0; i < domain().unitCell().nParameter(); ++i) {
         out << Int(i, 5)
             << "  "
             << Dbl(domain().unitCell().parameter(i), 18, 11)
             << std::endl;
      }
      out << std::endl;
   }

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
                                         domain().unitCell());
   }

   /*
   * Write w-fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(w_.hasData());
      domain_.fieldIo().writeFieldsRGrid(filename, w_.rgrid(),
                                         domain().unitCell(),
                                         w_.isSymmetric());
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
                                         domain().unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedGrid_);
      //UTIL_CHECK(hasCFields_);
      domain_.fieldIo().writeFieldsRGrid(filename, c_.rgrid(),
                                         domain().unitCell(),
                                         w_.isSymmetric());
   }

   /*
   * Write all concentration fields in real space (r-grid) format, for
   * each block (or solvent) individually rather than for each species.
   */
   template <int D>
   void System<D>::writeBlockCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(hasCFields_);

      // Create and allocate the DArray of fields to be written
      DArray< RField<D> > blockCFields;
      blockCFields.allocate(mixture_.nSolvent() + mixture_.nBlock());
      int n = blockCFields.capacity();
      for (int i = 0; i < n; i++) {
         blockCFields[i].allocate(domain().mesh().dimensions());
      }

      // Get data from Mixture and write to file
      mixture_.createBlockCRGrid(blockCFields);
      domain_.fieldIo().writeFieldsRGrid(filename, blockCFields,
                                         domain().unitCell(),
                                         w_.isSymmetric());
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
                                        domain().unitCell(),
                                        w_.isSymmetric());
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
                                        domain().unitCell(),
                                        w_.isSymmetric());
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
      domain().fieldIo().writeFieldHeader(file, 1, domain().unitCell(),
                                          w_.isSymmetric());
      file << "mesh " << std::endl
           << "          " << domain().mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
          file << "slice " << i << std::endl;
          domain().fieldIo().writeFieldRGrid(file, propagator.q(i),
                                             domain().unitCell(), 
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

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D>
   void System<D>::writeStars(std::string const & filename) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);
      bool isSymmetric = true;
      domain().fieldIo().writeFieldHeader(file, mixture_.nMonomer(),
                                          domain().unitCell(), 
                                          isSymmetric);
      domain_.basis().outputStars(file);
      file.close();
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D>
   void System<D>::writeWaves(std::string const & filename) const
   {
      UTIL_CHECK(domain_.hasGroup());
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);
      bool isSymmetric = true;
      domain().fieldIo().writeFieldHeader(file, mixture_.nMonomer(),
                                          domain().unitCell(), 
                                          isSymmetric);
      domain_.basis().outputWaves(file);
      file.close();
   }

   /*
   * Write all elements of the space group to a file.
   */
   template <int D>
   void System<D>::writeGroup(const std::string & filename) const
   {
      UTIL_CHECK(domain_.hasGroup());
      Pscf::Prdc::writeGroup(filename, domain_.group());
   }

   // Field format conversion functions

   /*
   * Convert fields from symmetry-adapted basis to real-space grid format.
   */
   template <int D>
   void System<D>::basisToRGrid(const std::string & inFileName,
                                const std::string & outFileName)
   {
      UTIL_CHECK(domain_.hasGroup());

      // If basis fields are not allocated, peek at field file header to
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
      }

      // Read, convert, and write fields
      UnitCell<D> tmpUnitCell;
      FieldIo<D> const & fieldIo = domain().fieldIo();
      fieldIo.readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo.convertBasisToRGrid(tmpFieldsBasis_, tmpFieldsRGrid_);
      fieldIo.writeFieldsRGrid(outFileName, tmpFieldsRGrid_,
                                 tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::rGridToBasis(const std::string & inFileName,
                                const std::string & outFileName)
   {
      UTIL_CHECK(domain_.hasGroup());

      // If basis fields are not allocated, peek at field file header to
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      FieldIo<D> const & fieldIo = domain().fieldIo();
      fieldIo.readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      fieldIo.convertRGridToBasis(tmpFieldsRGrid_, tmpFieldsBasis_);
      fieldIo.writeFieldsBasis(outFileName, tmpFieldsBasis_, 
                               tmpUnitCell);
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
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      FieldIo<D> const & fieldIo = domain().fieldIo();
      fieldIo.readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         domain().fft().inverseTransform(tmpFieldsKGrid_[i],
                                        tmpFieldsRGrid_[i]);
      }
      fieldIo.writeFieldsRGrid(outFileName, tmpFieldsRGrid_,
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
      if (domain_.hasGroup() && !isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      FieldIo<D> const & fieldIo = domain().fieldIo();
      fieldIo.readFieldsRGrid(inFileName, tmpFieldsRGrid_,
                                tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         domain().fft().forwardTransform(tmpFieldsRGrid_[i],
                                        tmpFieldsKGrid_[i]);
      }
      fieldIo.writeFieldsKGrid(outFileName, tmpFieldsKGrid_,
                                            tmpUnitCell);
   }

   /*
   * Convert fields from Fourier (k-grid) to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::kGridToBasis(const std::string & inFileName,
                                const std::string& outFileName)
   {
      UTIL_CHECK(domain_.hasGroup());

      // If basis fields are not allocated, peek at field file header to
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
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
      UTIL_CHECK(domain_.hasGroup());

      // If basis fields are not allocated, peek at field file header to
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
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
   * Check if r-grid fields have declared space group symmetry.
   */
   template <int D>
   bool
   System<D>::checkRGridFieldSymmetry(const std::string & inFileName,
                                      double epsilon)
   {
      UTIL_CHECK(domain_.hasGroup());

      // If basis fields are not allocated, peek at field file header to
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName);
         allocateFieldsBasis();
      }

      // Read fields
      UnitCell<D> tmpUnitCell;
      domain_.fieldIo().readFieldsRGrid(inFileName,
                                        tmpFieldsRGrid_, tmpUnitCell);

      // Check symmetry for all fields
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         bool symmetric;
         symmetric = domain_.fieldIo().hasSymmetry(tmpFieldsRGrid_[i],
                                                   epsilon);
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
      UTIL_CHECK(domain_.hasGroup());

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
      UTIL_CHECK(domain().mesh().size() > 0);
      UTIL_CHECK(!isAllocatedGrid_);

      // Alias for mesh dimensions
      IntVec<D> const & dimensions = domain().mesh().dimensions();

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
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(!isAllocatedBasis_);
      UTIL_CHECK(domain_.hasGroup());
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = domain_.basis().nBasis();
      UTIL_CHECK(nBasis > 0);

      w_.allocateBasis(nBasis);
      c_.allocateBasis(nBasis);

      // Temporary work space
      tmpFieldsBasis_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsBasis_[i].allocate(nBasis);
      }

      isAllocatedBasis_ = true;
   }

   /*
   * Peek at field file header, initialize unit cell parameters and basis.
   */
   template <int D>
   void System<D>::readFieldHeader(std::string filename)
   {
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(mixture_.nMonomer() > 0);

      // Open field file
      std::ifstream file;
      fileMaster_.openInputFile(filename, file);

      // Read field file header, and initialize basis if needed
      int nMonomer;
      bool isSymmetric;
      domain_.fieldIo().readFieldHeader(file, nMonomer,
                                        domain_.unitCell(), isSymmetric);
      // Function FieldIo::readFieldHeader initializes a basis if
      // hasGroup is true and the header contains a space group

      // Close field file
      file.close();

      // Postconditions
      UTIL_CHECK(mixture_.nMonomer() == nMonomer);
      UTIL_CHECK(domain().unitCell().nParameter() > 0);
      UTIL_CHECK(domain().unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(domain().unitCell().isInitialized());
      if (domain_.hasGroup() && isSymmetric) {
         UTIL_CHECK(domain_.basis().isInitialized());
         UTIL_CHECK(domain_.basis().nBasis() > 0);
      }

   }

   /*
   * Read a filename string and echo to log file (used in readCommands).
   */
   template <int D>
   void System<D>::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      if (in.fail()) {
          UTIL_THROW("Unable to read string parameter.");
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

      hasCFields_ = false;
      hasFreeEnergy_ = false;
   }

} // namespace Rpc
} // namespace Pscf
#endif
