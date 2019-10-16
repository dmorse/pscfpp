/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <pspg/GpuResources.h>

#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sys/time.h>
//#include <Windows.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <pspg/crystal/shiftToMinimum.h>
//#include <unistd.h>



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
      //basisPtr_(0),
      wavelistPtr_(0),
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
      hasFields_(0),
      hasSweep_(0)
   {  
      setClassName("System"); 
      std::cout<<"Constructor Called"<<std::endl;

      interactionPtr_ = new ChiInteraction();
      iteratorPtr_ = new FtsIterator<D>(this); 
      wavelistPtr_ = new WaveList<D>();
      //basisPtr_ = new Basis<D>();
      // sweepFactoryPtr_ = new SweepFactory(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   System<D>::~System()
   {
      delete interactionPtr_;
      delete iteratorPtr_; //there is an issue here. iterator ptr needs info of block size before initiation
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

      double time_start;
      double time_end;
      double time_wave_start;
      double time_wave_end;
      struct timeval tv;
      struct timezone tz;

      gettimeofday(&tv, &tz);
      time_start = (double) tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;

      std::cout<<"Read in mixture"<<std::endl;
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
      std::cout<<"Unit Cell read starting"<<std::endl;
      in >> unitCell_;
      hasUnitCell_ = true;
      std::cout<<"Unit Cell read completed"<<std::endl;

      IntVec<D> d;
      in >> d;
      mesh_.setDimensions(d);
      hasMesh_ = true;
      mixture().setMesh(mesh());
      mixture().setupUnitCell(unitCell());

      //std::string groupName;
      in >> groupName_;
      in >> groupName_;

      gettimeofday(&tv, &tz);
      time_wave_start = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;

      //basis().makeBasis(mesh(), unitCell(), groupName_);
      std::cout<<"Initializing wavelist"<<std::endl;
      wavelist().allocate(mesh(), unitCell());
      std::cout<<"Allocating completed"<<std::endl;
      wavelist().computeMinimumImages(mesh(), unitCell());
      std::cout<<"wavelist completed"<<std::endl;

      gettimeofday(&tv, &tz);
      time_wave_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      Log::file() << "wavelist initialized in  " 
                  << Dbl(time_wave_end - time_wave_start, 18, 11)<<'s' 
                  << std::endl;

      allocateFields();
      std::cout<<"fields allocated"<<std::endl;
      hasFields_ = true;


      // Initialize iterator
      readParamComposite(in, iterator());
      iterator().allocate();
      std::cout<<"Iterator Initialized"<<std::endl;

      gettimeofday(&tv, &tz);
      time_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      Log::file() << "System Parameters read in  " 
                  << Dbl(time_end - time_start, 18, 11)<<'s' 
                  << std::endl;
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
      double time_start;
      double time_end;
      struct timeval tv;
      struct timezone tz;

      gettimeofday(&tv, &tz);
      time_start = (double) tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;

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
         //wField(i).allocate(basis().nStar());
         wFieldGrid(i).allocate(mesh().dimensions());
         wFieldDft(i).allocate(mesh().dimensions());

         //cField(i).allocate(basis().nStar());
         cFieldGrid(i).allocate(mesh().dimensions());
         cFieldDft(i).allocate(mesh().dimensions());
      }

      //storageCFields_.allocate(mesh().dimensions());
      compositionKField_.allocate(mesh().dimensions());
      workArray.allocate(mesh().size());
   
      cudaMalloc((void**)&d_kernelWorkSpace_, NUMBER_OF_BLOCKS * sizeof(cufftReal));
      kernelWorkSpace_ = new cufftReal[NUMBER_OF_BLOCKS];
      //setup seed for rng
      timeval time;
      gettimeofday(&time, NULL);
      hasFields_ = true;

      gettimeofday(&tv, &tz);
      time_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      Log::file() << "System Containers allocated in  " 
                  << Dbl(time_end - time_start, 18, 11)<<'s' 
                  << std::endl;
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
         } 

         /*else if (command == "READ_WFIELDS") {
            in >> filename;
            Log::file() << " " << Str(filename, 20) <<std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(filename, inFile);
            readFields(inFile, wFields());
            inFile.close();

         } else if (command == "WRITE_WFIELDS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(filename, outFile);
            writeFields(outFile, wFields_);
            outFile.close();

            }*/ 
         else if (command == "WRITE_WFIELDGRIDS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(filename, outFile);
            //writeRFields(outFile, wFieldGrids_);
            writeRFields(outFile, wFieldGrids());
            outFile.close();

         }/* 
         else if (command == "WRITE_CFIELDS") {

            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(filename, outFile);
            writeFields(outFile, cFields_);
            outFile.close();

            }*/ 
         else if (command == "WRITE_CFIELDGRIDS") {
            in >> filename;
            Log::file() << "  " << Str(filename, 20) << std::endl;

            std::ofstream outFile;
            fileMaster().openOutputFile(filename, outFile);
            writeRFields(outFile, cFieldGrids_);
            outFile.close();
         } else if (command == "ITERATE") {
            Log::file() << std::endl;
            Log::file() << std::endl;

            double time_start;
            double time_end;
            struct timeval tv;
            struct timezone tz;

            gettimeofday(&tv, &tz);
            time_start = (double) tv.tv_sec + 
              (double)tv.tv_usec / 1000000.0;
            //input omega fields
            std::string inFileName;
            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);
            //changed this line
            readRFields(inFile, wFieldGrids());
            inFile.close();

            gettimeofday(&tv, &tz);
            time_end = (double)tv.tv_sec + 
               (double)tv.tv_usec / 1000000.0;
            Log::file() << "ReadFile_Time = " 
                        << Dbl(time_end - time_start, 18, 11)<<'s' 
                        << std::endl;

            //SYSTEMTIME timeStart, timeEnd;
            //double timeElapsed;
            //GetSystemTime(&timeStart);
            gettimeofday(&tv, &tz);
            time_start = (double) tv.tv_sec + 
              (double)tv.tv_usec / 1000000.0;
            std::cout<<"do we even reach here?"<<std::endl;
            int fail = iterator().solve();
            if(fail) {
               Log::file() << "Iterate has failed. Exiting "<<std::endl;
               exit(1);
            }
            //GetSystemTime(&timeEnd);
            //timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
            //timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
            //timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
            /*std::cout << " Time for SCFT ="
               << Dbl(timeElapsed, 18, 11) << 's' << std::endl;*/
            gettimeofday(&tv, &tz);
            time_end = (double)tv.tv_sec + 
               (double)tv.tv_usec / 1000000.0;
            computeFreeEnergy();
            outputThermo(Log::file());
            Log::file() << "SCF_Time = " 
            << Dbl(time_end - time_start, 18, 11)<<'s' 
            << std::endl;

         }
         /*else if (command == "FIELD_TO_RGRID") {
            std::string inFileName;
            std::string outFileName;

            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);
            readFields(inFile, cFields());
            inFile.close();

            //convert to rgrid
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               basis().convertFieldComponentsToDft(cField(i), cFieldDft(i));
               fft().inverseTransform(cFieldDft(i), cFieldGrid(i));
            }

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            writeRFields(outFile, cFieldGrids());
            outFile.close();

            } */
         /*else if (command == "RGRID_TO_FIELD") {
            std::string inFileName;
            std::string outFileName;

            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);
            readRFields(inFile, cFieldGrids());
            inFile.close();

            //convert to fields
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().forwardTransform(cFieldGrid(i), cFieldDft(i));
               basis().convertFieldDftToComponents(cFieldDft(i), cField(i));
            }

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            writeFields(outFile, cFields());
            outFile.close();

            }*/
         else if (command == "KGRID_TO_RGRID") {
            std::string inFileName;
            std::string outFileName;

            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) <<std::endl;

            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) <<std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);
            
            readKFields(inFile, cFieldDfts());
            inFile.close();

            //convert to rgrid
            for (int i = 0; i < mixture().nMonomer(); ++i) {
               fft().inverseTransform(cFieldDft(i), cFieldGrid(i));
            }

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            writeRFields(outFile, cFieldGrids());
            outFile.close();

         } else if (command == "RHO_TO_OMEGA") {

            std::string inFileName;
            std::string outFileName;

            in >> inFileName;
            Log::file() << " " << Str(inFileName, 20) << std::endl;

            in >> outFileName;
            Log::file() << " " << Str(outFileName, 20) << std::endl;

            std::ifstream inFile;
            fileMaster().openInputFile(inFileName, inFile);

            readRFields(inFile, cFieldGrids());
            
            inFile.close();

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

            std::ofstream outFile;
            fileMaster().openOutputFile(outFileName, outFile);
            writeRFields(outFile, wFieldGrids());
            outFile.close();

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
   template <int D>
   void System<D>::readFields(std::istream &in,
      DArray< RDField<D> >& fieldsToCopy)
   {
      UTIL_CHECK(hasMesh_);

      // Read grid dimensions
      std::string label;
      int nStar, nM;

      in >> label;
      UTIL_CHECK(label == "nStar");
      in >> nStar;
      UTIL_CHECK(nStar > 0);
      UTIL_CHECK(nStar == basis().nStar());

      in >> label;
      UTIL_CHECK(label == "nM");
      in >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == mixture().nMonomer());

      DArray<cufftReal*> fields;
      fields.allocate(nM);
      for (int i = 0; i < nM; i++) {
         fields[i] = new cufftReal[nStar];
      }
      // Read fields
      int i, j, idum;
      for (i = 0; i < nStar; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nM; ++j) {
            in >> fields[j][i];
         }
      }
      
      for (int i = 0; i < nM; i++) {
         cudaMemcpy(fieldsToCopy[i].cDField(), fields[i], sizeof(cufftReal) * nStar, cudaMemcpyHostToDevice);
         delete[] fields[i];
      }

   }
   */
   template <int D>
   void System<D>::readRFields(std::istream &in,
                                DArray<RDField<D> >& fields)
   {

      UTIL_CHECK(hasMesh_);

      std::string label;
      std::string groupName;
      IntVec<D> nGrid;
      int nM;
      int ver1, ver2;
      int dim;      
      FArray<double,6> params;

      in >> label;
      UTIL_CHECK(label == "format");
      in >> ver1;
      in >> ver2;
 
      in >> label;
      UTIL_CHECK(label == "dim");
      in >> dim;
      UTIL_CHECK(dim == D);

      unitCell().readHeader(in);
 
      in >> label;
      std::cout<<label<<std::endl;
      UTIL_CHECK(label == "group_name");
      in >> groupName;

      in >> label;
      std::cout<<label<<std::endl;
      UTIL_CHECK(label == "N_monomer");
      in >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == mixture().nMonomer());

      in >> label;
      UTIL_CHECK(label == "ngrid");
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());


      DArray<cufftReal*> temp;
      temp.allocate(nM);
      for(int i = 0; i < nM; ++i) {
         temp[i] = new cufftReal[mesh().size()];
      } 
      
      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         for(int k = 0; k < nM; ++k) {
            in >> std::setprecision(15)>> temp[k][rank];
         }
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         } 
      }
      
      for(int i = 0; i < nM; i++) {
         cudaMemcpy(fields[i].cDField(), temp[i],
            mesh().size() * sizeof(cufftReal), cudaMemcpyHostToDevice);
         delete[] temp[i];
         temp[i] = nullptr;
      }

   }

   template <int D>
   void System<D>::readKFields(std::istream &in,
                                 DArray<RDFieldDft<D> >& fields)
   {
      UTIL_CHECK(hasMesh_);

      std::string label;
      IntVec<D> nGrid;
      int nM;

      in >> label;
      UTIL_CHECK(label == "nGrid");
      in >> nGrid;
      UTIL_CHECK(nGrid == mesh().dimensions());

      in >> label;
      UTIL_CHECK(label == "nM");
      in >> nM;
      UTIL_CHECK(nM > 0);
      UTIL_CHECK(nM == mixture().nMonomer());

      // Read Fields;
      int kSize = 1;
      for (int i = 0; i < D; i++) {
        if (i == D - 1) {
           kSize *= (mesh().dimension(i) / 2 + 1);
        }
        else {
           kSize *= mesh().dimension(i);
        }
        
      }

      int idum;
      DArray<cufftComplex*> temp;
      temp.allocate(nM);
      for(int i = 0; i < nM; ++i) {
         temp[i] = new cufftComplex[kSize];
      }
      
      for(int i = 0; i < kSize; ++i) {
         in >> idum;
         for (int j = 0; j < nM; ++j) {
            in >> temp[j][i].x;
            in >> temp[j][i].y;
         }
      }
      
      for(int i = 0; i < nM; ++i) {
         cudaMemcpy(fields[i].cDField(), temp[i],
            kSize * sizeof(cufftComplex), cudaMemcpyHostToDevice);
         delete[] temp[i];
         temp[i] = nullptr;
      }
   }

   /*
   template <int D>
   void System<D>::writeFields(std::ostream &out,
      DArray<RDField<D> > const &  fieldsToCopy)
   {
      int i, j;
      int nStar = basis().nStar();
      int nM = mixture().nMonomer();
      out << "nStar     " << nStar << std::endl;
      out << "nM     " << nM << std::endl;

      DArray<cufftReal*> fields;
      fields.allocate(nM);
      for (int i = 0; i < nM; i++) {
         fields[i] = new cufftReal[nStar];
         cudaMemcpy(fields[i], fieldsToCopy[i].cDField(), nStar * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      }
      // Write fields
      for (i = 0; i < nStar; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nM; ++j) {
            out << "  " << Dbl(fields[j][i], 18, 11);
         }
         //out<< "  " << basis().wave(basis().star(i).beginId).indicesDft;
         out << std::endl;
      }
      for (int i = 0; i < nM; i++) {
         delete[] fields[i];
      }
   }
   */
   template <int D>
   void System<D>::writeRFields(std::ostream &out,
                           DArray<RDField<D> > const& fields)
   {
      double time_start;
      double time_end;
      struct timeval tv;
      struct timezone tz;

      gettimeofday(&tv, &tz);
      time_start = (double) tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;

      int nM = mixture().nMonomer();
      MeshIterator<D> itr(mesh().dimensions());
      //do not use white space like this...
      out << "format  1   0    " <<  std::endl;
      out << "dim    " <<  std::endl << "                    "<<D<< std::endl;
      unitCell().writeHeader(out);   
      out << "group_name    " <<  std::endl << "                    "<<groupName_<< std::endl;
      out << "N_monomer    " <<  std::endl << "                    "<<mixture().nMonomer()<< std::endl;
      out << "ngrid        " <<  std::endl<<"           "<<mesh().dimensions() << std::endl;

      DArray<cufftReal*> temp;
      temp.allocate(nM);
      for (int i = 0; i < nM; ++i) {
         temp[i] = new cufftReal[mesh().size()];
         cudaMemcpy(temp[i], fields[i].cDField(),
                    mesh().size() * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      }    

      IntVec<D> offsets;
      offsets[D - 1] = 1;
      for(int i = D - 1 ; i > 0; --i ) {
         offsets[i - 1] = offsets[i] * mesh().dimension(i);
      }
      IntVec<D> position;
      for(int i = 0; i < D; ++i) {
         position[i] = 0;
      }

      int rank = 0;
      int positionId;
      for(int i = 0; i < mesh().size(); i++) {
         rank = 0;
         for(int dim = 0; dim < D; ++dim) {
            rank += offsets[dim] * position[dim];
         }
         for(int k = 0; k < nM; ++k) {
            out << "  " << Dbl(temp[k][rank], 18, 15);
         }
         out<<'\n';
         //add position
         positionId = 0;
         while( positionId < D) {
            position[positionId]++;
            if ( position[positionId] == mesh().dimension(positionId) ) {
               position[positionId] = 0;
               positionId++;
               continue;
            }
            break;
         } 
      }
      
      for(int i = 0; i < nM; ++i) {
         delete[] temp[i];
         temp[i] = nullptr;
         
      }

      gettimeofday(&tv, &tz);
      time_end = (double)tv.tv_sec + 
         (double)tv.tv_usec / 1000000.0;
      Log::file() << "Files written in  " 
                  << Dbl(time_end - time_start, 18, 11)<<'s' 
                  << std::endl;
   }

   template <int D>
   void System<D>::writeKFields(std::ostream &out,
                           DArray<RDFieldDft<D> > const& fields)
   {
      int nM = mixture().nMonomer();
      MeshIterator<D> itr(mesh().dimensions());
      out << "nGrid   " << mesh().dimensions() << std::endl;
      out << "nM      " << nM                << std::endl;

      //write Fields
      DArray<cufftComplex*> temp;
     int kSize = 1;
     for (int i = 0; i < D; i++) {
        if (i == D - 1) {
           kSize *= (mesh().dimension(i) / 2 + 1);
        }
        else {
           kSize *= mesh().dimension(i);
        }
        
     }
      temp.allocate(nM);
      for(int i = 0; i < nM; ++i) {
         temp[i] = new cufftComplex[kSize];
         cudaMemcpy(temp[i], fields[i].cDField(), 
            kSize * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      }
     for (int i = 0; i < kSize; i++) {
        out << Int(i, 5);
        for (int j = 0; j < nM; ++j) {
           out << "  " << Dbl(temp[j][i].x, 18, 11)
              << Dbl(temp[j][i].y, 18, 11);
        }
        out << std::endl;
     }
      
      for(int i = 0; i < nM; ++i) {
         delete[] temp[i];
         temp[i] = nullptr;
      }
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
      //std::cout<<interaction().chi(0,1)<<std::endl;
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
