#ifndef PSPG_SYSTEM_TPP
#define PSPG_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <pspg/sweep/Sweep.h>
#include <pspg/sweep/SweepFactory.h>
#include <pspg/iterator/Iterator.h>
#include <pspg/iterator/IteratorFactory.h>
#include <pspg/field/RDField.h>
#include <pspg/math/GpuResources.h>

#include <pscf/inter/Interaction.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/homogeneous/Clump.h>

#include <util/param/BracketPolicy.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>
#include <getopt.h>

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
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      wavelistPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFieldsBasis_(),
      cFieldsBasis_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocated_(false),
      hasWFields_(false),
      hasCFields_(false),
      hasSymmetricFields_(false)
   {
      setClassName("System");
      domain_.setFileMaster(fileMaster_);

      interactionPtr_ = new Interaction();
      wavelistPtr_ = new WaveList<D>();
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this);
      sweepFactoryPtr_ = new SweepFactory<D>(*this);

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
      if (wavelistPtr_) {
         delete wavelistPtr_;
      }
      #if 0
      if (isAllocated_) {
         delete[] kernelWorkSpace_;
         cudaFree(d_kernelWorkSpace_);
      }
      #endif
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
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {

      readParamComposite(in, mixture());
      hasMixture_ = true;

      int nm = mixture().nMonomer();
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();

      // Initialize homogeneous object
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      // Read interaction (i.e., chi parameters)
      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain_);

      mixture().setMesh(mesh());

      // Construct wavelist
      wavelist().allocate(domain_.mesh(), unitCell());
      wavelist().computeMinimumImages(domain_.mesh(), unitCell());
      mixture().setupUnitCell(unitCell(), wavelist());

      // Allocate memory for w and c fields
      allocate();

      // Initialize iterator
      std::string className;
      bool isEnd;
      iteratorPtr_
            = iteratorFactoryPtr_->readObject(in, *this, className, isEnd);
      if (!iteratorPtr_) {
         std::string msg = "Unrecognized Iterator subclass name ";
         msg += className;
         UTIL_THROW(msg.c_str());
      }
      iterator().setup();

      // Optionally instantiate a Sweep object
      sweepPtr_ =
         sweepFactoryPtr_->readObjectOptional(in, *this, className, isEnd);
      if (sweepPtr_) {
         sweepPtr_->setSystem(*this);
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
   * Write parameter file, omitting the sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture_.writeParam(out);
      interaction().writeParam(out);
      domain_.writeParam(out);
      iteratorPtr_->writeParam(out);
      out << "}" << std::endl;
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
      UTIL_CHECK(isAllocated_);
      std::string command, filename, inFileName, outFileName;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command << std::endl;

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
         if (command == "COMPUTE") {
            compute();
         } else
         if (command == "ITERATE") {
            int error = iterate();
            if (error) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            sweep();
         } else
         if (command == "WRITE_W_BASIS") {
            UTIL_CHECK(hasWFields_);
            UTIL_CHECK(hasSymmetricFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsBasis(filename, wFieldsBasis(),
                                       unitCell());
         } else
         if (command == "WRITE_W_RGRID") {
            UTIL_CHECK(hasWFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, wFieldsRGrid(),
                                       unitCell());
         } else
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(hasCFields_);
            UTIL_CHECK(hasSymmetricFields_);
            fieldIo().writeFieldsBasis(filename, cFieldsBasis(), 
                                       unitCell());
         } else
         if (command == "WRITE_C_RGRID") {
            UTIL_CHECK(hasCFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, cFieldsRGrid(), 
                                       unitCell());
         } else
         if (command == "WRITE_C_BLOCK_RGRID") {
            readEcho(in, filename);
            writeBlockCRGrid(filename);
         } else
         if (command == "WRITE_PROPAGATOR_TAIL") {
            readEcho(in, filename);
            int polymerID, blockID;
            in >> polymerID;
            in >> blockID;
            Log::file() << Str("polymer ID   ", 21) << polymerID << "\n"
                        << Str("block ID   ", 21) << blockID << std::endl;
            writePropagatorTail(filename, polymerID, blockID);
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
            // Read from file in k-grid format
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToRGrid(inFileName, outFileName);
         } else
         if (command == "GUESS_W_FROM_C") {
            // Read c field file in r-grid format
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            guessWfromC(inFileName, outFileName);
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
      if (np > 0) {
         for (i = 0; i < np; ++i) {

            // Initial array of clump sizes
            for (j = 0; j < nm; ++j) {
               c_[j] = 0.0;
            }

            // Compute clump sizes for all monomer types.
            nb = mixture_.polymer(i).nBlock();
            for (k = 0; k < nb; ++k) {
               Block<D>& block = mixture_.polymer(i).block(k);
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

   /*
   * Compute Helmholtz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(hasWFields_);
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
             (wFieldsRGrid_[i].cDField(), cFieldsRGrid_[i].cDField(),
              workArray.cDField(), nx);
         temp += gpuSum(workArray.cDField(),nx) / double(nx);
      }
      fHelmholtz_ -= temp;

      // Compute excess interaction free energy
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
           assignUniformReal<<<nBlocks, nThreads>>>
               (workArray.cDField(), interaction().chi(i, j), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>
               (workArray.cDField(), cFieldsRGrid_[i].cDField(), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>
               (workArray.cDField(), cFieldsRGrid_[j].cDField(), nx);
           fHelmholtz_ += gpuSum(workArray.cDField(), nx) / double(nx);
         }
      }

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


   template <int D>
   void System<D>::writeThermo(std::ostream& out)
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
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

   // W-Field Modifier Functions

   /*
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      fieldIo().readFieldsBasis(filename, wFieldsBasis_,
                                domain_.unitCell());
      fieldIo().convertBasisToRGrid(wFieldsBasis(), wFieldsRGrid_);
      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;
   }

   /*
   * Read w-field in real space grid (r-grid) format.
   */
   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      // Update system wFields
      int nMonomer = mixture().nMonomer();
      int nBasis = basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> & w = wFieldsBasis_[i];
         for (int j = 0; j < nBasis; ++j) {
            w[j] = f[j];
         }
      }

      // Update system wFieldsRgrid
      domain_.fieldIo().convertBasisToRGrid(wFieldsBasis_, wFieldsRGrid_);

      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray<Field> const & fields)
   {
      int nMonomer = mixture_.nMonomer();
      int nMesh = domain_.mesh().size();
      for (int i = 0; i < nMonomer; ++i) {
         wFieldsRGrid_[i] = fields[i];
      }
      hasWFields_ = true;
      hasSymmetricFields_ = false;
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using array of r-grid fields as input.
   */
   template <int D>
   void System<D>::setWRGrid(DField<cudaReal> & fields)
   {
      const int nMonomer = mixture_.nMonomer();
      const int nMesh = domain_.mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      //cudaReal const * src;
      //cudaReal* dst;
      for (int i = 0; i < nMonomer; i++) {
         //dst = wFieldsRGrid_[i].cDField();
         //src = fields.cDField() + i*nMesh;
         //assignReal<<<nBlocks, nThreads>>>(dst, src, nMesh);
         assignReal<<<nBlocks, nThreads>>>(wFieldsRGrid_[i].cDField(), 
                                           fields.cDField() + i*nMesh, 
                                           nMesh);
         //cudaMemcpy(dst, src, nMesh*sizeof(cudaReal),
         //           cudaMemcpyDeviceToDevice);
      }

      hasWFields_ = true;
      hasSymmetricFields_ = false;
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using array of r-grid fields as input.
   */
   template <int D>
   void System<D>::symmetrizeWFields()
   {
      UTIL_CHECK(hasWFields_);
      fieldIo().convertRGridToBasis(wFieldsRGrid(), wFieldsBasis_);
      fieldIo().convertBasisToRGrid(wFieldsBasis(), wFieldsRGrid_);
      hasSymmetricFields_ = true;
      hasCFields_ = false;
   }

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      UTIL_CHECK(domain_.unitCell().lattice() == unitCell.lattice());
      domain_.unitCell() = unitCell;
      mixture_.setupUnitCell(unitCell, wavelist());
      wavelist().computedKSq(domain_.unitCell());
   }

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      UTIL_CHECK(domain_.unitCell().nParameter() == parameters.size());
      domain_.unitCell().setParameters(parameters);
      mixture_.setupUnitCell(domain_.unitCell(), wavelist());
      wavelist().computedKSq(domain_.unitCell());
   }

   // Primary SCFT Computations

   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {
      fieldIo().readFieldsRGrid(filename, wFieldsRGrid_,
                                 domain_.unitCell());
      hasWFields_ = true;
      hasSymmetricFields_ = false;
      hasCFields_ = false;
   }

   /*
   * Solve MDE for current w-fields, without iteration.
   */
   template <int D>
   void System<D>::compute(bool needStress)
   {
      UTIL_CHECK(hasWFields_);

      // Solve the modified diffusion equation (without iteration)
      mixture().compute(wFieldsRGrid(), cFieldsRGrid_);
      hasCFields_ = true;

      // Convert c fields from r-grid to basis format
      if (hasSymmetricFields_) {
         fieldIo().convertRGridToBasis(cFieldsRGrid(), cFieldsBasis_);
      }

      if (needStress) {
         mixture().computeStress(wavelist());
      }
   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   template <int D>
   int System<D>::iterate()
   {
      UTIL_CHECK(hasWFields_);
      hasCFields_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator
      int error = iterator().solve();

      hasCFields_ = true;

      if (hasSymmetricFields_) {
         fieldIo().convertRGridToBasis(cFieldsRGrid(), cFieldsBasis_);
      }

      if (!error) {
         if (!iterator().isFlexible()) {
            mixture().computeStress(wavelist());
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
      UTIL_CHECK(hasWFields_);
      UTIL_CHECK(hasSymmetricFields_);
      UTIL_CHECK(sweepPtr_);
      UTIL_CHECK(hasSweep());
      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
   }

   /*
   * Write w-fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeWBasis(const std::string & filename)
   {
      UTIL_CHECK(hasWFields_);
      UTIL_CHECK(hasSymmetricFields_);
      fieldIo().writeFieldsBasis(filename, wFieldsBasis(), unitCell());
   }

   /*
   * Write w-fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasWFields_);
      fieldIo().writeFieldsRGrid(filename, wFieldsRGrid(), unitCell());
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(const std::string & filename)
   {
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(hasSymmetricFields_);
      fieldIo().writeFieldsBasis(filename, cFieldsBasis(), unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);
      fieldIo().writeFieldsRGrid(filename, cFieldsRGrid_, unitCell());
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
      DArray<Field> blockCFields;
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
   * Write the last time slice of the propagator.
   */
   template <int D>
   void System<D>::writePropagatorTail(const std::string & filename,
                                       int polymerID, int blockID)
   {
      const cudaReal* d_tailField = 
                      mixture_.polymer(polymerID).propagator(blockID, 1).tail();

      RDField<D> tailField;
      tailField.allocate(mesh().size());
      cudaMemcpy(tailField.cDField(), d_tailField,
                 mesh().size() * sizeof(cudaReal), cudaMemcpyDeviceToDevice);
      fieldIo().writeFieldRGrid(filename, tailField, unitCell());
   }

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D>
   void System<D>::writeStars(const std::string & outFileName) const
   {
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(),
                                 unitCell());
      basis().outputStars(outFile);
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D>
   void System<D>::writeWaves(const std::string & outFileName) const
   {
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(), 
                                 unitCell());
      basis().outputWaves(outFile);
   }

   // Field Operations

   /*
   * Convert fields from symmetry-adpated basis to real-space grid format.
   */
   template <int D>
   void System<D>::basisToRGrid(const std::string & inFileName,
                                const std::string & outFileName)
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToRGrid(tmpFieldsBasis_, tmpFieldsRGrid_);
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_, tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::rGridToBasis(const std::string & inFileName,
                                const std::string & outFileName)
   {
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
                                const std::string& outFileName) const
   {
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
                                const std::string & outFileName) const
   {
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
                                const std::string& outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      fieldIo().convertKGridToBasis(tmpFieldsKGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert fields from symmetry-adapted basis to Fourier (k-grid) format.
   */
   template <int D>
   void System<D>::basisToKGrid(const std::string & inFileName,
                                const std::string & outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToKGrid(tmpFieldsBasis_, tmpFieldsKGrid_);
      fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, tmpUnitCell);
   }

   #if 0
   /*
   * Check if fields are symmetric under space group.
   */
   template <int D>
   bool System<D>::checkRGridFieldSymmetry(const std::string & inFileName) 
   const
   {
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

   /*
   * Construct guess for omega (w-field) from rho (c-field).
   *
   * Modifies wFields and wFieldsRGrid and also outputs wFields to file.
   */
   template <int D>
   void System<D>::guessWfromC(std::string const & inFileName, 
                               std::string const & outFileName)
   {
      // Read c fields from input file
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, 
                                domain_.unitCell());

      // Compute w fields from c fields
      for (int i = 0; i < basis().nBasis(); ++i) {
         for (int j = 0; j < mixture_.nMonomer(); ++j) {
            wFieldsBasis_[j][i] = 0.0;
            for (int k = 0; k < mixture_.nMonomer(); ++k) {
               wFieldsBasis_[j][i] += interaction().chi(j,k) 
                                      * tmpFieldsBasis_[k][i];
            }
         }
      }

      // Convert to r-grid format
      fieldIo().convertBasisToRGrid(wFieldsBasis_, wFieldsRGrid_);
      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;

      // Write w field in basis format
      fieldIo().writeFieldsBasis(outFileName, wFieldsBasis(), unitCell());

      #if 0
      // Older code from G.K. Cheong's version
      fieldIo().readFieldsRGrid(inFileName, cFieldsRGrid_,
                                domain_.unitCell());

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(mesh().size(), nBlocks, nThreads);

      // Compute w fields, excluding Lagrange multiplier contribution
      //code is bad here, `mangled' access of data in array
      for (int i = 0; i < mixture().nMonomer(); ++i) {
         assignUniformReal<<<nBlocks, nThreads>>>
             (wFieldsRGrid_[i].cDField(), 0, mesh().size());
      }
      for (int i = 0; i < mixture().nMonomer(); ++i) {
         for (int j = 0; j < mixture().nMonomer(); ++j) {
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (wFieldsRGrid_[i].cDField(), cFieldsRGrid_[j].cDField(),
                 interaction().chi(i,j), mesh().size());
         }
      }

      // Write w fields to file in r-grid format
      fieldIo().writeFieldsRGrid(outFileName, wFieldsRGrid(),
                                 unitCell());
      #endif

   }

   // Private functions

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocate()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFieldsBasis_.allocate(nMonomer);
      wFieldsRGrid_.allocate(nMonomer);

      cFieldsBasis_.allocate(nMonomer);
      cFieldsRGrid_.allocate(nMonomer);

      tmpFieldsBasis_.allocate(nMonomer);
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);

      int nBasis = basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         wFieldsBasis_[i].allocate(nBasis);
         wFieldsRGrid_[i].allocate(mesh().dimensions());

         cFieldsBasis_[i].allocate(nBasis);
         cFieldsRGrid_[i].allocate(mesh().dimensions());

         tmpFieldsBasis_[i].allocate(nBasis);
         tmpFieldsRGrid_[i].allocate(mesh().dimensions());
         tmpFieldsKGrid_[i].allocate(mesh().dimensions());
      }

      workArray.allocate(mesh().size());
      ThreadGrid::setThreadsLogical(mesh().size());

      isAllocated_ = true;
   }

} // namespace Pspg
} // namespace Pscf
#endif
