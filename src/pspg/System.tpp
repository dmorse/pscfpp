#ifndef PSPG_SYSTEM_TPP
#define PSPG_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <pspg/math/GpuResources.h>

#include <pscf/homogeneous/Clump.h>
#include <pscf/crystal/shiftToMinimum.h>

#include <pspg/iterator/IteratorMediatorCUDA.h>
#include <pspg/iterator/IteratorFactory.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

//#include <iomanip>
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
      iteratorMediatorPtr_(0),
      iteratorFactoryPtr_(0),
      wavelistPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocated_(false),
      hasWFields_(false),
      hasCFields_(false)
      // hasSweep_(0)
   {  
      setClassName("System");
      domain_.setFileMaster(fileMaster_);

      interactionPtr_ = new ChiInteraction();
      iteratorMediatorPtr_ = new IteratorMediatorCUDA<D>(*this);
      iteratorFactoryPtr_ = new IteratorFactory<D>(*iteratorMediatorPtr_); 
      wavelistPtr_ = new WaveList<D>();

      ThreadGrid::init();
      
      // sweepFactoryPtr_ = new SweepFactory(*this);
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
      if (iteratorMediatorPtr_) {
         delete iteratorMediatorPtr_;
      }
      if (iteratorFactoryPtr_) {
         delete iteratorFactoryPtr_;
      }
      if (wavelistPtr_) {
         delete wavelistPtr_; 
      }
      if (isAllocated_) {
         delete[] kernelWorkSpace_;
         cudaFree(d_kernelWorkSpace_);
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
      bool tFlag = false;  // GPU input threads (maximum number of threads per block)
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
            iFlag = true;
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
      //int ns = mixture().nSolvent(); 
      int ns = 0;

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
      wavelist().allocate(mesh(), unitCell());
      wavelist().computeMinimumImages(mesh(), unitCell());
      mixture().setupUnitCell(unitCell(), wavelist());

      // Allocate memory for w and c fields
      allocate();

      // Initialize iterator through the factory and mediator
      std::string className;
      bool isEnd;
      Iterator<FieldCUDA>* iteratorPtr
         = iteratorFactoryPtr_->readObject(in, *this, className, isEnd);
      if (!iteratorPtr) {
         std::string msg = "Unrecognized Iterator subclass name ";
         msg += className;
         UTIL_THROW(msg.c_str());
      }
      iteratorMediator().setIterator(*iteratorPtr);
      iteratorMediator().setup();
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
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocate()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      wFieldsRGrid_.allocate(nMonomer);
      wFieldsKGrid_.allocate(nMonomer);

      cFields_.allocate(nMonomer);
      cFieldsRGrid_.allocate(nMonomer);
      cFieldsKGrid_.allocate(nMonomer);

      tmpFields_.allocate(nMonomer);
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);

      for (int i = 0; i < nMonomer; ++i) {
         wField(i).allocate(basis().nStar());
         wFieldRGrid(i).allocate(mesh().dimensions());
         wFieldKGrid(i).allocate(mesh().dimensions());

         cField(i).allocate(basis().nStar());
         cFieldRGrid(i).allocate(mesh().dimensions());
         cFieldKGrid(i).allocate(mesh().dimensions());

         tmpFields_[i].allocate(basis().nStar());
         tmpFieldsRGrid_[i].allocate(mesh().dimensions());
         tmpFieldsKGrid_[i].allocate(mesh().dimensions());
      }

      //storageCFields_.allocate(mesh().dimensions());
      //compositionKField_.allocate(mesh().dimensions());
      workArray.allocate(mesh().size());
      ThreadGrid::setThreadsLogical(mesh().size());
   
      cudaMalloc((void**)&d_kernelWorkSpace_, ThreadGrid::nBlocks() * sizeof(cudaReal));
      kernelWorkSpace_ = new cudaReal[ThreadGrid::nBlocks()];

      isAllocated_ = true;
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

            fieldIo().readFieldsBasis(filename, wFields());
            fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
            hasWFields_ = true;

         } else
         if (command == "COMPUTE") {
            // Read w (chemical potential fields) if not done previously 
            if (!hasWFields_) {
               readEcho(in, filename);
               readWBasis(filename);
            }
            // Solve the modified diffusion equation, without iteration
            compute();
         } else 
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);

            fieldIo().readFieldsRGrid(filename, wFieldsRGrid());
            hasWFields_ = true;

         } else 
         if (command == "ITERATE") {
            // Read w (chemical potential) fields if not done previously 
            if (!hasWFields_) {
               readEcho(in, filename);
               readWBasis(filename);
            }
            // Attempt iteration to convergence
            int fail = iterate();
            if (fail) {
               readNext = false;
            }
         } else 
         if (command == "WRITE_W_BASIS") {
            UTIL_CHECK(hasWFields_);
            readEcho(in, filename);
            fieldIo().convertRGridToBasis(wFieldsRGrid(), wFields());
            fieldIo().writeFieldsBasis(filename, wFields());
         } else 
         if (command == "WRITE_W_RGRID") {
            UTIL_CHECK(hasWFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, wFieldsRGrid());
         } else 
         if (command == "WRITE_C_BASIS") {
            UTIL_CHECK(hasCFields_);
            readEcho(in, filename);
            fieldIo().convertRGridToBasis(cFieldsRGrid(), cFields());
            fieldIo().writeFieldsBasis(filename, cFields());
         } else 
         if (command == "WRITE_C_RGRID") {
            UTIL_CHECK(hasCFields_);
            readEcho(in, filename);
            fieldIo().writeFieldsRGrid(filename, cFieldsRGrid());
         } else 
         if (command == "WRITE_PROPAGATOR") {
            int polymerID, blockID;
            readEcho(in, filename);
            in >> polymerID;
            in >> blockID;
            Log::file() << Str("polymer ID   ", 21) << polymerID << "\n"
                        << Str("block ID   ", 21) << blockID << std::endl;
            writePropagatorRGrid(filename, polymerID, blockID);
         } else
         if (command == "BASIS_TO_RGRID") {
            hasCFields_ = false;
            readEcho(in, inFileName);
            readEcho(in, outFileName);

            fieldIo().readFieldsBasis(inFileName, cFields());
            fieldIo().convertBasisToRGrid(cFields(), cFieldsRGrid());
            fieldIo().writeFieldsRGrid(outFileName, cFieldsRGrid());
         } else 
         if (command == "RGRID_TO_BASIS") {
            hasCFields_ = false;
            readEcho(in, inFileName);
            readEcho(in, outFileName);

            fieldIo().readFieldsRGrid(inFileName, cFieldsRGrid());
            fieldIo().convertRGridToBasis(cFieldsRGrid(), cFields());
            fieldIo().writeFieldsBasis(outFileName, cFields());

         } else 
         if (command == "KGRID_TO_RGRID") {
            hasCFields_ = false;

            // Read from file in k-grid format
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            fieldIo().readFieldsKGrid(inFileName, cFieldsKGrid());

            // Use FFT to convert k-grid r-grid
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().inverseTransform(cFieldKGrid(i), cFieldRGrid(i));
            }
            // Write to file in r-grid format
            fieldIo().writeFieldsRGrid(outFileName, cFieldsRGrid());

         } else 
         if (command == "RHO_TO_OMEGA") {

            // GPU resources
            int nBlocks, nThreads;
            ThreadGrid::setThreadsLogical(mesh().size(), nBlocks, nThreads);

            // Read c field file in r-grid format
            readEcho(in, inFileName);
            readEcho(in, outFileName);

            fieldIo().readFieldsRGrid(inFileName, cFieldsRGrid());

            // Compute w fields, excluding Lagrange multiplier contribution
            //code is bad here, `mangled' access of data in array
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               assignUniformReal<<<nBlocks, nThreads>>>(wFieldRGrid(i).cDField(), 0, mesh().size());
            }
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               for (int j = 0; j < mixture().nMonomer(); ++j) {
                  pointWiseAddScale<<<nBlocks, nThreads>>>(wFieldRGrid(i).cDField(), cFieldRGrid(j).cDField(), 
                                                            interaction().chi(i,j), mesh().size());
               }
            }

            // Write w fields to file in r-grid format
            fieldIo().writeFieldsRGrid(outFileName, wFieldsRGrid());

         } else {

            Log::file() << "  Error: Unknown command  " << command << std::endl;
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
            Block<D>& block = mixture().polymer(i).block(k);
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

   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      fHelmholtz_ = 0.0;
 
      // Compute ideal gas contributions to fHelhmoltz_
      Polymer<D>* polymerPtr;
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

      int nm  = mixture().nMonomer();
      int nx = mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nx, nBlocks, nThreads);
      
      double temp = 0;
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
           assignUniformReal<<<nBlocks, nThreads>>>(workArray.cDField(), interaction().chi(i, j), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>(workArray.cDField(), cFieldsRGrid_[i].cDField(), nx);
           inPlacePointwiseMul<<<nBlocks, nThreads>>>(workArray.cDField(), cFieldsRGrid_[j].cDField(), nx);
           fHelmholtz_ += (gpuSum(workArray.cDField(), nx) / nx);
         }
         
         assignReal<<<nBlocks, nThreads>>>(workArray.cDField(), wFieldsRGrid_[i].cDField(), nx);
         inPlacePointwiseMul<<<nBlocks, nThreads>>>(workArray.cDField(), cFieldsRGrid_[i].cDField(), nx);
         temp += gpuSum(workArray.cDField(), nx);
      }
      fHelmholtz_ -= (temp / nx);
      
      // Compute pressure
      pressure_ = -fHelmholtz_;
      for (int i = 0; i < np; ++i) {
         polymerPtr = &mixture().polymer(i);
         phi = polymerPtr->phi();
         mu = polymerPtr->mu();
         length = polymerPtr->length();
         pressure_ += mu * phi /length;
      }

   }


   template <int D>
   void System<D>::outputThermo(std::ostream& out)
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

   /*  
   * Read w-field in symmetry adapted basis format.
   */  
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {   
      fieldIo().readFieldsBasis(filename, wFields());
      fieldIo().convertBasisToRGrid(wFields(), wFieldsRGrid());
      domain_.basis().update();
      hasWFields_ = true;
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
      mixture().compute(wFieldsRGrid(), cFieldsRGrid());

      // Convert c and w fields from r-grid to basis
      fieldIo().convertRGridToBasis(wFieldsRGrid(), wFields());
      fieldIo().convertRGridToBasis(cFieldsRGrid_, cFields_);

      hasCFields_ = true;

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
      int error = iteratorMediator().solve();
      
      hasCFields_ = true;
      fieldIo().convertRGridToBasis(wFieldsRGrid(), wFields());
      fieldIo().convertRGridToBasis(cFieldsRGrid(), cFields());

      if (!error) {   
         if (!domain().isFlexible()) {
            mixture().computeStress(wavelist());
         }
         computeFreeEnergy();
         outputThermo(Log::file());
      }
      return error;
   }

   /*
   * Convert fields from symmetry-adpated basis to real-space grid format.
   */
   template <int D>
   void System<D>::basisToRGrid(const std::string & inFileName,
                                const std::string & outFileName)
   {
      fieldIo().readFieldsBasis(inFileName, tmpFields_);
      fieldIo().convertBasisToRGrid(tmpFields_, tmpFieldsRGrid_);
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_);
   }

   /*
   * Write w-fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeWBasis(const std::string & filename)
   {
      UTIL_CHECK(hasWFields_);
      fieldIo().writeFieldsBasis(filename, wFields());
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(const std::string & filename)
   {
      UTIL_CHECK(hasCFields_);
      fieldIo().writeFieldsBasis(filename, cFields());
   }

   template <int D>
   void System<D>::writePropagatorRGrid(const std::string & filename, int polymerID, int blockID)
   {
      const cudaReal* d_tailField = mixture_.polymer(polymerID).propagator(blockID, 1).tail();
      
      // convert this cudaReal pointer to an RDField. Yikes.
      RDField<D> tailField;
      tailField.allocate(mesh().size());
      cudaMemcpy(tailField.cDField(), d_tailField, mesh().size() * sizeof(cudaReal), cudaMemcpyDeviceToDevice);
      // output. 
      fieldIo().writeFieldRGrid(filename, tailField);
   }

   /*  
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */  
   template <int D>
   void System<D>::rGridToBasis(const std::string & inFileName,
                                const std::string & outFileName)
   {   
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_);
      fieldIo().convertRGridToBasis(tmpFieldsRGrid_, tmpFields_);
      fieldIo().writeFieldsBasis(outFileName, tmpFields_);
   }

} // namespace Pspg
} // namespace Pscf
#endif 
