#ifndef PSPG_SYSTEM_TPP
#define PSPG_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <pspg/GpuResources.h>

#include <pscf/homogeneous/Clump.h>
#include <pscf/crystal/shiftToMinimum.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>
#include <getopt.h>

//#include <unistd.h>
//#include <Windows.h>

//global variable for kernels
int THREADS_PER_BLOCK;
int NUMBER_OF_BLOCKS;

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
      mesh_(),
      unitCell_(),
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      basisPtr_(),
      //fft_(),
      //groupName_(),
      wavelistPtr_(0),
      fieldIo_(),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(0),
      hasUnitCell_(0),
      hasFields_(0)
      // hasSweep_(0)
   {  
      setClassName("System"); 

      interactionPtr_ = new ChiInteraction();
      iteratorPtr_ = new AmIterator<D>(this); 
      wavelistPtr_ = new WaveList<D>();
      basisPtr_ = new Basis<D>();

      // sweepFactoryPtr_ = new SweepFactory(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   System<D>::~System()
   {
      delete interactionPtr_;
      delete iteratorPtr_; 
      delete wavelistPtr_; 
      delete[] kernelWorkSpace_;
      cudaFree(d_kernelWorkSpace_);
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
      bool wFlag = false;  // GPU input 1 (# of blocks)
      bool tFlag = false;  // GPU input 2 (threads per block)
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f1:2:")) != -1) {
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
         case '1': //number of blocks
            NUMBER_OF_BLOCKS = atoi(optarg);
            wFlag = true;
            break;
         case '2': //threads per block
            THREADS_PER_BLOCK = atoi(optarg);
            tFlag = true;
            //something like this
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

      if (!wFlag) {
         std::cout<<"Number of blocks not set " <<std::endl;
         exit(1);
      }

      if (!tFlag) {
         std::cout<<"Threads per block not set " <<std::endl;
         exit(1);
      }

   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {

      #if 0
      double time_start;
      double time_end;
      double time_wave_start;
      double time_wave_end;
      struct timeval tv;
      struct timezone tz;

      gettimeofday(&tv, &tz);
      time_start = (double) tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      #endif

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

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());
      // ----------------------

      read(in, "unitCell", unitCell_);
      hasUnitCell_ = true;
      
      // IntVec<D> d;
      read(in, "mesh", mesh_);
      hasMesh_ = true;
      mixture().setMesh(mesh());

      read(in, "groupName", groupName_);

      #if 0
      in >> d;
      mesh_.setDimensions(d);
      hasMesh_ = true;
      mixture().setMesh(mesh());
      #endif


      #if 0
      gettimeofday(&tv, &tz);
      time_wave_start = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      #endif

      // Construct wavelist 
      wavelist().allocate(mesh(), unitCell());
      wavelist().computeMinimumImages(mesh(), unitCell());
      mixture().setupUnitCell(unitCell(), wavelist());

      basis().makeBasis(mesh(), unitCell(), groupName_);

      //std::cout<<"nstar = "<<basis().nStar();

      fieldIo_.associate(unitCell_, mesh_, fft_, groupName_,
                         *basisPtr_, fileMaster_);

      #if 0
      gettimeofday(&tv, &tz);
      time_wave_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      Log::file() << "wavelist initialized in  " 
                 << Dbl(time_wave_end - time_wave_start, 18, 11) << 's' 
                 << std::endl;
      #endif

      allocateFields();
      hasFields_ = true;

      // Initialize iterator
      readParamComposite(in, iterator());
      iterator().allocate();

      //std::cout<<"nstar = "<<basis().nStar();

      #if 0
      gettimeofday(&tv, &tz);
      time_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
         Log::file() << "System Parameters read in  " 
                  << Dbl(time_end - time_start, 18, 11)<<'s' 
                  << std::endl;
      #endif
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
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::allocateFields()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasMesh_);

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      wFieldGrids_.allocate(nMonomer);
      wFieldDfts_.allocate(nMonomer);

      cFields_.allocate(nMonomer);
      cFieldGrids_.allocate(nMonomer);
      cFieldDfts_.allocate(nMonomer);

      //size of grid is based on basis function
      for (int i = 0; i < nMonomer; ++i) {
         wField(i).allocate(basis().nStar());
         wFieldGrid(i).allocate(mesh().dimensions());
         wFieldDft(i).allocate(mesh().dimensions());

         cField(i).allocate(basis().nStar());
         cFieldGrid(i).allocate(mesh().dimensions());
         cFieldDft(i).allocate(mesh().dimensions());
      }

      //storageCFields_.allocate(mesh().dimensions());
      compositionKField_.allocate(mesh().dimensions());
      workArray.allocate(mesh().size());
   
      cudaMalloc((void**)&d_kernelWorkSpace_, NUMBER_OF_BLOCKS * sizeof(cufftReal));
      kernelWorkSpace_ = new cufftReal[NUMBER_OF_BLOCKS];
      hasFields_ = true;

   }

   
   /*
   * Read and execute commands from a specified command file.
   */
   //will add more commands as they are tested
   template <int D>
   void System<D>::readCommands(std::istream &in) 
   {
      UTIL_CHECK(hasFields_);
      std::string command;
      std::string filename;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command <<std::endl;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else 
         if (command == "READ_WFIELDS") {
            in >> filename;
            Log::file() << " " << Str(filename, 20) <<std::endl;

            fieldIo().readFieldsBasis(filename, wFields());

         } else 
         if (command == "WRITE_W_BASIS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsBasis(filename, wFields());
         } else 
         if (command == "WRITE_WFIELDGRIDS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsRGrid(filename, wFieldGrids());
         } else 
         if (command == "WRITE_C_BASIS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsBasis(filename, cFields());
         } 
         else if (command == "WRITE_CFIELDGRIDS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;
            fieldIo().writeFieldsRGrid(filename, cFieldGrids_);
         } 
         else if (command == "ITERATE") {
            Log::file() << std::endl;
            Log::file() << std::endl;

            // Read w fields
            in >> filename;
            Log::file() << " " << Str(filename, 20) <<std::endl;
            fieldIo().readFieldsRGrid(filename, wFieldGrids());

            int fail = iterator().solve();

            if (!fail) {
               computeFreeEnergy();
               outputThermo(Log::file());
            } else {
               Log::file() << "Iterate has failed. Exiting "<<std::endl;
               exit(1);
            }
         } else 
         if (command == "BASIS_TO_RGRID") {
            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;


            fieldIo().readFieldsBasis(inFileName, cFields());
            fieldIo().convertBasisToRGrid(cFields(), cFieldGrids());

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsRGrid(outFileName, cFieldGrids());
         } else 
         if (command == "RGRID_TO_BASIS") {
            std::string inFileName;
            std::string outFileName;

            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;

            fieldIo().readFieldsRGrid(inFileName, cFieldGrids());
            fieldIo().convertRGridToBasis(cFieldGrids(), cFields());
            fieldIo().writeFieldsBasis(outFileName, cFields());

         } else 
         if (command == "KGRID_TO_RGRID") {

            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;
            fieldIo().readFieldsKGrid(inFileName, cFieldDfts());

            // Convert to rgrid
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().inverseTransform(cFieldDft(i), cFieldGrid(i));
            }

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;
            fieldIo().writeFieldsRGrid(outFileName, cFieldGrids());

         } else 
         if (command == "RHO_TO_OMEGA") {

            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) << std::endl;
            fieldIo().readFieldsRGrid(inFileName, cFieldGrids());

            //code is bad here, `mangled' access of data in array
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               assignUniformReal << < NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (wFieldGrid(i).cDField(), 0, mesh().size());
            }
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               for (int j = 0; j < mixture().nMonomer(); ++j) {
                  pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (wFieldGrid(i).cDField(), cFieldGrid(j).cDField(), 
                     interaction().chi(i,j), mesh().size());
               }
            }

            std::string outFileName;
            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) << std::endl;
            fieldIo().writeFieldsRGrid(outFileName, wFieldGrids());

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
         std::cout<< mu <<std::endl;
         std::cout << length<<std::endl;
         fHelmholtz_ += phi*( mu - 1.0 )/length;
      }

      int nm  = mixture().nMonomer();
      int nx = mesh().size();
      //RDField<D> workArray;
      //workArray.allocate(nx);
      float temp = 0;
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
           assignUniformReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(workArray.cDField(), interaction().chi(i, j), nx);
           inPlacePointwiseMul << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (workArray.cDField(), cFieldGrids_[i].cDField(), nx);
           inPlacePointwiseMul << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (workArray.cDField(), cFieldGrids_[j].cDField(), nx);
           fHelmholtz_ += (reductionH(workArray, nx) / nx);
         }
         
         assignReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (workArray.cDField(), wFieldGrids_[i].cDField(), nx);
         inPlacePointwiseMul << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (workArray.cDField(), cFieldGrids_[i].cDField(), nx);
         temp += reductionH(workArray, nx);
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

      //out << "Polymers:" << std::endl;
      //out << "    i"
      //    << "        phi[i]      "
      //    << "        mu[i]       " 
      //    << std::endl;
      for (int i = 0; i < mixture().nPolymer(); ++i) {
         out << Int(i, 5) 
             << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
             << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)  
             << std::endl;
      }
      out << std::endl;
   }

   template <int D>
   cufftReal System<D>::innerProduct(const RDField<D>& a, const RDField<D>& b, int size) {

     switch(THREADS_PER_BLOCK){
     case 512:
       deviceInnerProduct<512><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 256:
       deviceInnerProduct<256><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 128:
       deviceInnerProduct<128><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 64:
       deviceInnerProduct<64><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 32:
       deviceInnerProduct<32><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 16:
       deviceInnerProduct<16><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 8:
       deviceInnerProduct<8><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 4:
       deviceInnerProduct<4><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 2:
       deviceInnerProduct<2><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     case 1:
       deviceInnerProduct<1><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_kernelWorkSpace_, a.cDField(), b.cDField(), size);
       break;
     }
      cudaMemcpy(kernelWorkSpace_, d_kernelWorkSpace_, NUMBER_OF_BLOCKS * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      cufftReal final = 0;
      cufftReal c = 0;
      //use kahan summation to reduce error
      for (int i = 0; i < NUMBER_OF_BLOCKS; ++i) {
         cufftReal y = kernelWorkSpace_[i] - c;
         cufftReal t = final + y;
         c = (t - final) - y;
         final = t;

      }

      return final;
   }

   template<int D>
   cufftReal System<D>::reductionH(const RDField<D>& a, int size) {
     reduction << < NUMBER_OF_BLOCKS / 2, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal) >> > (d_kernelWorkSpace_, a.cDField(), size);
      cudaMemcpy(kernelWorkSpace_, d_kernelWorkSpace_, NUMBER_OF_BLOCKS / 2 * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      cufftReal final = 0;
      cufftReal c = 0;
      for (int i = 0; i < NUMBER_OF_BLOCKS / 2; ++i) {
         cufftReal y = kernelWorkSpace_[i] - c;
         cufftReal t = final + y;
         c = (t - final) - y;
         final = t;
         }
      return final;
   }

} // namespace Pspg
} // namespace Pscf
#endif
