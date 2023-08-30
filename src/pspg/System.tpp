#ifndef PSPG_SYSTEM_TPP
#define PSPG_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <pspg/compressor/Compressor.h>
#include <pspg/compressor/CompressorFactory.h>
#include <pspg/sweep/Sweep.h>
#include <pspg/sweep/SweepFactory.h>
#include <pspg/iterator/Iterator.h>
#include <pspg/iterator/IteratorFactory.h>
#include <pspg/field/RDField.h>
#include <pspg/math/GpuResources.h>

#include <pscf/inter/Interaction.h>
#include <pscf/math/IntVec.h>
#include <pscf/homogeneous/Clump.h>

#include <util/param/BracketPolicy.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

#include <string>
#include <unistd.h>
//#include <getopt.h>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   System<D>::System()
    : mixture_(),
      domain_(),
      mcSimulator_(*this),
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
      fHelmholtz_(0.0),
      fIdeal_(0.0),
      fInter_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      hasMcSimulator_(false),
      isAllocatedGrid_(false),
      isAllocatedBasis_(false),
      hasCFields_(false)
   {
      setClassName("System");
      domain_.setFileMaster(fileMaster_);
      w_.setFieldIo(domain_.fieldIo());

      interactionPtr_ = new Interaction();
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this);
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      compressorFactoryPtr_ = new CompressorFactory<D>(*this);
      BracketPolicy::set(BracketPolicy::Optional);
      ThreadGrid::init();
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
      bool eflag = false; // echo
      bool pFlag = false; // param file
      bool cFlag = false; // command file
      bool iFlag = false; // input prefix
      bool oFlag = false; // output prefix
      bool tFlag = false; // GPU input threads (max. # of threads per block)
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
      int tArg = 0;

      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:t:")) != -1) {
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
         case 't': //number of threads per block, user set
            tArg = atoi(optarg);
            tFlag = true;
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

      // If option -t, set the threads per block.
      if (tFlag) {
         ThreadGrid::setThreadsPerBlock(tArg);
      }

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
   {  readParam(fileMaster().paramFile()); }
   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {

      // Read the Mixture{ ... } block
      readParamComposite(in, mixture());
      hasMixture_ = true;

      int nm = mixture().nMonomer();
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(np >= 0);
      UTIL_CHECK(ns >= 0);
      UTIL_CHECK(np + ns > 0);

      // Read the Interaction{ ... } block
      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      // Read the Domain{ ... } block
      readParamComposite(in, domain_);

      mixture().setMesh(mesh());
      mixture().setupUnitCell(unitCell());
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(domain_.unitCell().nParameter() > 0);
      UTIL_CHECK(domain_.unitCell().lattice() != UnitCell<D>::Null);

      // Allocate memory for w and c fields
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
         Log::file() << indent() << "  [Iterator{} absent]\n";
      }

      // Optionally instantiate a Sweep object (if an iterator exists)
      if (iteratorPtr_) {
         sweepPtr_ = 
            sweepFactoryPtr_->readObjectOptional(in, *this, className, 
                                                 isEnd);
         if (!sweepPtr_ && ParamComponent::echo()) {
            Log::file() << indent() << "  Sweep{ [absent] }\n";
         }
      }
      
      // Optionally instantiate a Compressor object
      compressorPtr_ = 
         compressorFactoryPtr_->readObjectOptional(in, *this, className, 
                                                   isEnd);
      if (!compressorPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Compressor{ [absent] }\n";
      }
      
      // Optionally read an McSimulator
      readParamCompositeOptional(in, mcSimulator_);
      if (mcSimulator_.isActive()) {
         hasMcSimulator_ = true;
      } else {
         hasMcSimulator_ = false;
      }
      

      // Initialize homogeneous object
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

   }

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
            // Read c field file in r-grid format
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
            // Solve the modified diffusion equation for current w fields
            compute();
         } else
         if (command == "ITERATE") {
            // Attempt to iteratively solve a SCFT problem
            int error = iterate();
            if (error) {
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
            UTIL_CHECK(hasCompressor());
            compressor().compress();
         } else
         if (command == "SIMULATE") {
            // Perform a field theoretic MC simulation
            int nStep;
            in >> nStep;
            Log::file() << "   "  << nStep << "\n";
            simulate(nStep);
         } else
         if (command == "ANALYZE_TRAJECTORY") {
            int min;
            in >> min;
            int max;
            in >> max;
            std::string classname;
            readEcho(in, classname);
            readEcho(in, filename);
            mcSimulator_.analyzeTrajectory(min, max, classname, filename);
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
         if (command == "WRITE_W_BASIS") {
            UTIL_CHECK(w_.hasData());
            UTIL_CHECK(w_.isSymmetric());
            readEcho(in, filename);
            fieldIo().writeFieldsBasis(filename, w_.basis(),
                                       domain_.unitCell());
         } else
         if (command == "WRITE_W_RGRID") {
            UTIL_CHECK(w_.hasData());
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, w_.rgrid(),
                                       domain_.unitCell());
         } else
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(hasCFields_);
            UTIL_CHECK(w_.isSymmetric());
            fieldIo().writeFieldsBasis(filename, c_.basis(), 
                                       domain_.unitCell());
         } else
         if (command == "WRITE_C_RGRID") {
            UTIL_CHECK(hasCFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, c_.rgrid(), 
                                       domain_.unitCell());
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
                    << Str("segment ID  ", 21) << segmentId << std::endl;
            writeQSlice(filename, polymerId, blockId, directionId, 
                                  segmentId);
         } else
         if (command == "WRITE_Q_TAIL") {
            int polymerId, blockId, directionId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                      << Str("block ID  ", 21) << blockId << "\n"
                      << Str("direction ID  ", 21) << directionId << "\n";
            writeQTail(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q") {
            int polymerId, blockId, directionId;
            readEcho(in, filename);
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
            //writeGroup(filename, domain_.group());
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
         if (command == "KGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToBasis(inFileName, outFileName);
         } else
         if (command == "BASIS_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            basisToKGrid(inFileName, outFileName);
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
   template <int D>
   void System<D>::readCommands()
   {
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile());
   }

   // W-Field Modifier Functions

   /*
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(filename); 
         allocateFieldsBasis(); 
      }

      // Read w fields
      w_.readBasis(filename, domain_.unitCell());
      domain_.waveList().computeKSq(domain_.unitCell());
      domain_.waveList().computedKSq(domain_.unitCell());
      mixture_.setupUnitCell(domain_.unitCell(), domain_.waveList());
      hasCFields_ = false;
      hasMcHamiltonian_ = false;
   }

   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(filename);  
         allocateFieldsBasis();  
      }

      w_.readRGrid(filename, domain_.unitCell());
      domain_.waveList().computeKSq(domain_.unitCell());
      domain_.waveList().computedKSq(domain_.unitCell());
      mixture_.setupUnitCell(domain_.unitCell(), domain_.waveList());
      hasCFields_ = false;
      hasMcHamiltonian_ = false;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      w_.setBasis(fields);
      hasCFields_ = false;
      hasMcHamiltonian_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray< RDField<D> > const & fields)
   {
      w_.setRGrid(fields);
      hasCFields_ = false;
      hasMcHamiltonian_ = false;
   }

   /*
   * Set new w-field values, using unfoldeded array of r-grid fields.
   */
   template <int D>
   void System<D>::setWRGrid(DField<cudaReal> & fields)
   {
      w_.setRGrid(fields);
      hasCFields_ = false;
   }

   /*
   * Construct estimate for w-fields from c-fields.
   *
   * Modifies wFields and wFieldsRGrid.
   */
   template <int D>
   void System<D>::estimateWfromC(std::string const & filename)
   {
      UTIL_CHECK(hasMixture_);
      const int nMonomer = mixture_.nMonomer();
      UTIL_CHECK(nMonomer > 0);

      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(filename); 
         allocateFieldsBasis(); 
      }
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = domain_.basis().nBasis();
      UTIL_CHECK(nBasis > 0);

      // Allocate temporary storage
      DArray< DArray<double> > tmpCFieldsBasis;
      tmpCFieldsBasis.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpCFieldsBasis[i].allocate(nBasis);
      }

      // Read c fields from input file
      fieldIo().readFieldsBasis(filename, tmpCFieldsBasis, 
                                domain_.unitCell());

      // Compute w fields from c fields
      for (int i = 0; i < nBasis; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            tmpFieldsBasis_[j][i] = 0.0;
            for (int k = 0; k < nMonomer; ++k) {
               tmpFieldsBasis_[j][i] += interaction().chi(j,k) 
                                        * tmpCFieldsBasis[k][i];
            }
         }
      }

      // Store estimated w fields in System w container
      w_.setBasis(tmpFieldsBasis_);

      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using array of r-grid fields as input.
   */
   template <int D>
   void System<D>::symmetrizeWFields()
   {  w_.symmetrize(); }

   // Unit Cell Modifiers

   /*
   * Set the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      domain_.setUnitCell(unitCell);
      mixture_.setupUnitCell(unitCell, domain_.waveList());
      if (!isAllocatedBasis_) {
         allocateFieldsBasis();
      }
   }

   /*
   * Set the associated unit cell.
   */
   template <int D>
   void 
   System<D>::setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                          FSArray<double, 6> const & parameters)
   {
      domain_.setUnitCell(lattice, parameters);
      mixture_.setupUnitCell(domain_.unitCell(), domain_.waveList());
      if (!isAllocatedBasis_) {
         allocateFieldsBasis();
      }
   }

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      domain_.setUnitCell(parameters);
      mixture_.setupUnitCell(domain_.unitCell(), domain_.waveList());
      if (!isAllocatedBasis_) {
         allocateFieldsBasis();
      }
   }

   // Primary SCFT Computations

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
      mixture().compute(w_.rgrid(), c_.rgrid());
      hasCFields_ = true;
      // Convert c fields from r-grid to basis format
      if (w_.isSymmetric()) {
         fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis());
      }

      if (needStress) {
         mixture().computeStress(domain_.waveList());
      }
   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   template <int D>
   int System<D>::iterate(bool isContinuation)
   {
      UTIL_CHECK(w_.hasData());
      hasCFields_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator
      int error = iterator().solve(isContinuation);

      hasCFields_ = true;

      if (w_.isSymmetric()) {
         fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis());
      }

      if (!error) {
         if (!iterator().isFlexible()) {
            mixture().computeStress(domain_.waveList());
         }
         computeFreeEnergy();
         writeThermo(Log::file());
      }
      hasMcHamiltonian_ = false;
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
      UTIL_CHECK(sweepPtr_);
      UTIL_CHECK(hasSweep());
      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
   }
   
   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void System<D>::simulate(int nStep)
   {
      UTIL_CHECK(hasCompressor());
      UTIL_CHECK(hasMcSimulator_);
      mcSimulator_.simulate(nStep);
   }

   // Thermodynamic Properties

   /*
   * Compute Helmholtz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(hasCFields_);

      // Initialize to zero
      fHelmholtz_ = 0.0;

      double phi, mu;
      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      // Compute polymer ideal gas contributions to fHelhmoltz_
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double phi, mu, length;
         int np = mixture().nPolymer();
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

      // Compute solvent ideal gas contributions to fHelhmoltz_
      if (ns > 0) {
         Solvent<D>* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture_.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1E-08) {
               fHelmholtz_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      int nm  = mixture().nMonomer();
      int nx = mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nx, nBlocks, nThreads);

      // Compute Legendre transform subtraction
      double temp = 0.0;
      for (int i = 0; i < nm; i++) {
         pointWiseBinaryMultiply<<<nBlocks,nThreads>>>
             (w_.rgrid(i).cDField(), c_.rgrid()[i].cDField(),
              workArray_.cDField(), nx);
         temp += gpuSum(workArray_.cDField(),nx) / double(nx);
      }
      fHelmholtz_ -= temp;
      fIdeal_ = fHelmholtz_;

      // Compute excess interaction free energy
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
           assignUniformReal<<<nBlocks, nThreads>>>
               (workArray_.cDField(), interaction().chi(i, j), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>
               (workArray_.cDField(), c_.rgrid()[i].cDField(), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>
               (workArray_.cDField(), c_.rgrid()[j].cDField(), nx);
           fHelmholtz_ += gpuSum(workArray_.cDField(), nx) / double(nx);
         }
      }
      fInter_ = fHelmholtz_ - fIdeal_;

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

   // Output Operations

   /*
   * Write parameter file, omitting the sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture_.writeParam(out);
      interaction().writeParam(out);
      domain_.writeParam(out);
      if (iteratorPtr_) {
         iteratorPtr_->writeParam(out);
      }
      out << "}" << std::endl;
   }

   template <int D>
   void System<D>::writeThermo(std::ostream& out)
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
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
      for (int i = 0; i < unitCell().nParameter(); ++i) {
         out << Int(i, 5)
             << "  " 
             << Dbl(unitCell().parameter(i), 18, 11)
             << std::endl;
      }
      out << std::endl;
   }

   /*
   * Write w-fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeWBasis(const std::string & filename)
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      fieldIo().writeFieldsBasis(filename, w_.basis(), unitCell());
   }

   /*
   * Write w-fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(const std::string & filename) const
   {
      UTIL_CHECK(w_.hasData());
      fieldIo().writeFieldsRGrid(filename, w_.rgrid(), unitCell());
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(const std::string & filename)
   {
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(w_.isSymmetric());
      fieldIo().writeFieldsBasis(filename, c_.basis(), unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);
      fieldIo().writeFieldsRGrid(filename, c_.rgrid(), unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format, for each
   * block (or solvent) individually rather than for each species.
   */
   template <int D>
   void System<D>::writeBlockCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);

      // Create and allocate the DArray of fields to be written
      DArray< RDField<D> > blockCFields;
      blockCFields.allocate(mixture_.nSolvent() + mixture_.nBlock());
      int n = blockCFields.capacity();
      for (int i = 0; i < n; i++) {
         blockCFields[i].allocate(mesh().dimensions());
      }

      // Get data from Mixture and write to file
      mixture_.createBlockCRGrid(blockCFields);
      fieldIo().writeFieldsRGrid(filename, blockCFields, unitCell());
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
      Propagator<D> const& 
          propagator = polymer.propagator(blockId, directionId);
      RDField<D> field;
      field.allocate(mesh().size());
      cudaMemcpy(field.cDField(), propagator.q(segmentId),
               mesh().size() * sizeof(cudaReal), cudaMemcpyDeviceToDevice);
      fieldIo().writeFieldRGrid(filename, field, unitCell());
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
      Propagator<D> const& 
          propagator = polymer.propagator(blockId, directionId);
      RDField<D> field;
      field.allocate(mesh().size());
      cudaMemcpy(field.cDField(), propagator.tail(),
                 mesh().size()*sizeof(cudaReal), 
                 cudaMemcpyDeviceToDevice);
      fieldIo().writeFieldRGrid(filename, field, unitCell());
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
      Propagator<D> const& 
           propagator = polymer.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);

      // Write header
      fieldIo().writeFieldHeader(file, 1, unitCell());
      file << "ngrid" << std::endl
           << "          " << mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      RDField<D> field;
      field.allocate(mesh().size());
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
          file << "slice " << i << std::endl;
          cudaMemcpy(field.cDField(), propagator.q(i),
                     mesh().size() * sizeof(cudaReal), 
                     cudaMemcpyDeviceToDevice);
          fieldIo().writeFieldRGrid(file, field, unitCell(), hasHeader);
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
   void System<D>::writeStars(const std::string & filename) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream outFile;
      fileMaster_.openOutputFile(filename, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(),
                                 unitCell());
      basis().outputStars(outFile);
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D>
   void System<D>::writeWaves(const std::string & filename) const
   {
      UTIL_CHECK(domain_.basis().isInitialized());
      std::ofstream outFile;
      fileMaster_.openOutputFile(filename, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(), 
                                 unitCell());
      basis().outputWaves(outFile);
   }

   /*
   * Write all elements of the space group to a file.
   */
   template <int D>
   void System<D>::writeGroup(const std::string & filename) const
   {  
      Pscf::writeGroup(filename, domain_.group()); 
   }

   // Field File Operations

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
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
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
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      fieldIo().convertRGridToBasis(tmpFieldsRGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, 
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
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         fft().inverseTransform(tmpFieldsKGrid_[i], tmpFieldsRGrid_[i]);
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
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, 
                                tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         fft().forwardTransform(tmpFieldsRGrid_[i], tmpFieldsKGrid_[i]);
      }
      fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, 
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
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      fieldIo().convertKGridToBasis(tmpFieldsKGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, 
                                 tmpUnitCell);
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
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }

      // Read, convert and write fields
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToKGrid(tmpFieldsBasis_, tmpFieldsKGrid_);
      fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, 
                                 tmpUnitCell);
   }

   // Private member functions

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocateFieldsGrid()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      int nMonomer = mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(!isAllocatedGrid_);

      // Alias for mesh dimensions
      IntVec<D> const & dimensions = domain_.mesh().dimensions();

      // Allocate W Fields
      w_.setNMonomer(nMonomer);
      w_.allocateRGrid(dimensions);

      // Allocate C Fields
      c_.setNMonomer(nMonomer);
      c_.allocateRGrid(dimensions);

      // Allocate temporary work space
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsRGrid_[i].allocate(dimensions);
         tmpFieldsKGrid_[i].allocate(dimensions);
      }

      workArray_.allocate(mesh().size());
      ThreadGrid::setThreadsLogical(mesh().size());

      isAllocatedGrid_ = true;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocateFieldsBasis()
   {
      // Preconditions and constants
      UTIL_CHECK(hasMixture_);
      const int nMonomer = mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(isAllocatedGrid_);
      UTIL_CHECK(!isAllocatedBasis_);
      UTIL_CHECK(domain_.unitCell().nParameter() > 0);
      UTIL_CHECK(domain_.mesh().size() > 0);
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = basis().nBasis();
      UTIL_CHECK(nBasis > 0);

      w_.allocateBasis(nBasis);
      c_.allocateBasis(nBasis);

      // Temporary work space
      tmpFieldsBasis_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsBasis_[i].allocate(nBasis);
      }

      // Allocate memory for waveList
      if (!domain_.waveList().hasMinimumImages()) {
         domain_.waveList().computeMinimumImages(domain_.mesh(), 
                                                 domain_.unitCell());
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
      domain_.fieldIo().readFieldHeader(file, nMonomer, 
                                        domain_.unitCell());
      // FieldIo::readFieldHeader initializes a basis if needed
      file.close();

      // Postconditions
      UTIL_CHECK(mixture_.nMonomer() == nMonomer);
      UTIL_CHECK(domain_.unitCell().nParameter() > 0);
      UTIL_CHECK(domain_.unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(domain_.unitCell().isInitialized());
      UTIL_CHECK(domain_.basis().isInitialized());
      UTIL_CHECK(domain_.basis().nBasis() > 0);
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
   * Initialize Pscf::Homogeneous::Mixture homogeneous_ member.
   */
   template <int D>
   void System<D>::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer();
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();

      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      int i;   // molecule index
      int j;   // monomer index

      // Loop over polymer molecule species
      if (np > 0) {

         DArray<double> cTmp;
         cTmp.allocate(nm);

         int k;   // block or clump index
         int nb;  // number of blocks
         int nc;  // number of clumps
         for (i = 0; i < np; ++i) {

            // Initial array of clump sizes
            for (j = 0; j < nm; ++j) {
               cTmp[j] = 0.0;
            }

            // Compute clump sizes for all monomer types.
            nb = mixture_.polymer(i).nBlock();
            for (k = 0; k < nb; ++k) {
               Block<D>& block = mixture_.polymer(i).block(k);
               j = block.monomerId();
               cTmp[j] += block.length();
            }

            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (cTmp[j] > 1.0E-8) {
                  ++nc;
               }
            }
            homogeneous_.molecule(i).setNClump(nc);

            // Set clump properties for this Homogeneous::Molecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (cTmp[j] > 1.0E-8) {
                  homogeneous_.molecule(i).clump(k).setMonomerId(j);
                  homogeneous_.molecule(i).clump(k).setSize(cTmp[j]);
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
      hasMcHamiltonian_ = false;
   }

   #if 0
   /*
   * Check if fields are symmetric under space group.
   */
   template <int D>
   bool System<D>::checkRGridFieldSymmetry(const std::string & inFileName) 
   const
   {
      // If basis fields are not allocated, peek at field file header to 
      // get unit cell parameters, initialize basis and allocate fields.
      if (!isAllocatedBasis_) {
         readFieldHeader(inFileName); 
         allocateFieldsBasis(); 
      }
      UTIL_CHECK(domain_.basis().isInitialized());
      const int nBasis = domain_.basis().nBasis();
      UTIL_CHECK(nBasis > 0);

      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         bool symmetric = fieldIo().hasSymmetry(tmpFieldsRGrid_[i]);
         if (!symmetric) {
            return false;
         }
      }
      return true;
   }
   #endif

} // namespace Pspg
} // namespace Pscf
#endif
