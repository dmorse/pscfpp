#ifndef PSSP_GPU_BASIS_TPP
#define PSSP_GPU_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "Basis.h"
#include "cuComplex.h"
#include <pssp_gpu/GpuResources.h>

static __global__ void makeDksqHelper(cufftReal* dksq, const int* waveBz,
                                      const cufftReal* dkkBasis,
                                      int nParams, int nStar, int dim) {
   //actual size is nStar*nParams
   //each thread do two calculation
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   //I do not think this is a good order but I just copied the cpu code for now
   for(int param = 0; param < nParams; ++param) {
      for (int i = startID; i < nStar; i += nThreads) {
         dksq[(param * nStar) + i] = 0;
      }
   }
        
   for(int param = 0; param < nParams; ++param) {
      for (int i = startID; i < nStar; i += nThreads) {
          for(int j = 0; j < dim; ++j) {
            for(int k = 0; k < dim; ++k) {
               dksq[(param * nStar) + i] += waveBz[i*dim + j] 
                  * waveBz[i*dim +k]
                  * dkkBasis[k + (j * dim) + (param * dim * dim)];
            }
         }
      }
   }
}

namespace Pscf {
namespace Pssp_gpu
{
   template <int D>
   Basis<D>::Basis()
      : nWave_(0), nStar_(0), unitCellPtr_(0), mesh_(0)
   {}

   template <int D>
   void Basis<D>::makeBasis(const Mesh<D>& mesh, const UnitCell<D>& unitCell,
                            std::string groupName)
   {
      if (groupName == "I") {
         //To Do: Need some method to obtain the space group symmetry 
         //       if other than I
      } else {
         UTIL_THROW("Unimplemented space group");
      }

      //Associate both mesh and unit cell
      mesh_ = &mesh;
      unitCellPtr_ = &unitCell;
      
      MeshIterator<D> itr(mesh.dimensions());
      
      //make waves
      waves_.allocate(mesh.size());
      waveId_.allocate(mesh.size());
      nWave_ = mesh.size();

      for (itr.begin(); !itr.atEnd(); ++itr) {
         waves_[itr.rank()].sqNorm = unitCell.ksq(itr.position());
         waves_[itr.rank()].indicesDft = itr.position();

         if ((itr.position(D-1) + 1) > (mesh_->dimension(D-1)/2 + 1)) {
            waves_[itr.rank()].implicit = true;
         } else {
            waves_[itr.rank()].implicit = false;
         }

         //`unsorted' waves; waves appear in grid order
         waveId_[itr.rank()] = itr.rank();
      }
      
      //make stars
      //To do: If not I, sort according to increasing ksq values
      //Changed: The size of stars_ is not known apriori. Used a GArray
      //For `I', reserved the right size.
      stars_.reserve(mesh.size());
      waveBz_ = new int[mesh.size() * D];
      cudaMalloc((void**)&waveBz_d, sizeof(int) * mesh.size() * D);

      int beginWave = 0;
      int endWave = 0;
      int iStar = 0;
      bool cancel;
      int invertFlag = 0;
      std::complex<double> cNorm;

      for (itr.begin(); !itr.atEnd() ; ++itr) {
         stars_[itr.rank()].invertFlag = 3;
      }

      for (itr.begin(); !itr.atEnd() ; ++itr) {
         //To do: If not I, create a list of waves with equal ksq 
         //magnitudes. The I space group implicitly assumes each waves belong
         //in a single star and the sorting is not done.
         
         endWave = beginWave+1;
         double starPhase[1] = {0};
         //To do: unimplemented local scope function
         //cancel = isCancelled(waves_[beginWave].indices);
         cancel = false;

         #if 0
         //check for cancellation consistency in star
         for(int j = beginWave; j < endWave; ++j) {
            if( cancel != isCancelled(waves_[j].indices))
               UTIL_THROW("Inconsistent cancellation in star");
         }
         #endif

         if (!cancel) {
            if (invertFlag == 1) { 
               cNorm = exp(std::complex<double>(0,1) *
                          starPhase[endWave-1-beginWave]);
            } else { 
               cNorm = exp(std::complex<double>(0,1)*starPhase[0]);
            }

            cNorm *= sqrt(double(endWave - beginWave));
         }

         for (int j = beginWave; j < endWave; ++j) {
            waves_[j].starId = iStar;
            if (!cancel) {
               waves_[j].coeff = exp(std::complex<double>(0,1) *
                                    starPhase[j - beginWave])/cNorm;
               //std::cout<<waves_[j].coeff<<std::endl;
            } else {   
               waves_[j].coeff = exp(std::complex<double>(0,0));
            }
         }
         
         //fill up stars_ object
         //unimplemented: sign flag.
         stars_[iStar].size = endWave - beginWave;
         stars_[iStar].beginId = beginWave;
         stars_[iStar].endId = endWave - 1;
         stars_[iStar].cancel = cancel;
         ++nStar_;

         //To do: a method here to decide what is the actual invertFlag
         //A lot easier when waves are sorted in min(ksq). Currently, the second
         //wave of 001 index 1 is paired with say 009.
         IntVec<D> G1 = itr.position();
         IntVec<D> G2;

         bool isClosed = true;
         for (int j = 0; j < D; ++j) {
            G2[j] = -G1[j];
         }
         mesh_->shift(G2);
         if (G2 != G1) {
            isClosed = false;
         }
         if (isClosed) {
            invertFlag = 0;
            stars_[iStar].invertFlag = invertFlag;
         } else {
            mesh_->shift(G2);
            int partnerId = mesh_->rank(G2);
            //std::cout<<stars_[partnerId].invertFlag<<std::endl;
            if (stars_[partnerId].invertFlag == 3) {
               stars_[iStar].invertFlag = 1;
            } else {
               stars_[iStar].invertFlag = -1;
            }
         }

         /*
         if(itr.position(D-1) > (mesh_.dimensions(D-1)/2 + 1) )
            stars_[iStar].invertFlag = -1;
         else
            stars_[iStar].invertFlag = 1;
         */

         if (invertFlag == -1) {
            shiftToMinimum(waves_[endWave-1].indicesDft, mesh_->dimensions(), waveBz_ + (iStar * D));
         } else {
            shiftToMinimum(waves_[beginWave].indicesDft, mesh_->dimensions(), waveBz_ + (iStar * D));
         }

         iStar++;
         beginWave = endWave;
      }

      cudaMemcpy(waveBz_d, waveBz_, sizeof(int) * mesh.size() * D, cudaMemcpyHostToDevice);
      //make lookup table for gpu conversions
      IntVec<D> offsets;
      IntVec<D> G2;
      int partnerId;
      offsets[D - 1] = 1;
      for (int j = D - 1; j > 0; --j) {
         if (j == D - 1) {
            offsets[j - 1] = offsets[j] * (mesh_->dimension(j) / 2 + 1);
         }
         else {
            offsets[j - 1] = offsets[j] * mesh_->dimension(j);
         }
      }

      //allocate space for table
      int* rankTable = new int[mesh_->size()];
      int* partnerIdTable = new int[mesh_->size()];
      bool* implicitTable = new bool[mesh_->size()];
      cufftComplex* coeffTable = new cufftComplex[mesh_->size()];
      int* starIdTable = new int[mesh_->size()];
      int* invertFlagTable = new int[mesh_->size()];
      int* star2rank = new int[nStar_];

      //allocate space on device side
	  //std::cout << "monkey out" << std::endl;
      cudaMalloc((void**) &rankTable_, sizeof(int) * mesh_->size());
      cudaMalloc((void**) &partnerIdTable_, sizeof(int) * mesh_->size());
      cudaMalloc((void**) &implicitTable_, sizeof(bool) * mesh_->size());
      cudaMalloc((void**) &coeffTable_, sizeof(cufftComplex) * mesh_->size());
      cudaMalloc((void**) &starIdTable_, sizeof(int) * mesh_->size());
	  //std::cout << "what is wrong?" << std::endl;
      cudaMalloc((void**) &invertFlagTable_, sizeof(int) * mesh_->size());
      cudaMalloc((void**) &star2rank_, sizeof(int) * nStar_);
	  //leaking memory need to delete somewhere
      cudaMalloc((void**) &starIdTable_, sizeof(int) * mesh_->size());


      for (int i = 0; i < mesh_->size(); ++i) {
         IntVec<D> waveId = mesh_->position(i);

         //fill implicitTable
         implicitTable[i] = wave(waveId).implicit;

         //fill partnerIdTable
         for(int j = 0; j < D; ++j) {
            G2[j] = -waveId[j];
         }
         mesh_->shift(G2);
         partnerId = mesh_->rank(G2);
         partnerIdTable[i] = partnerId;

         //fill rankTable
         int rank = 0;
         for(int j = 0; j < D; ++j) {
            rank += wave(waveId).indicesDft[j] * offsets[j];
         }
         rankTable[i] = rank;

         //fill coeff table
         coeffTable[i].x = wave(waveId).coeff.real();
         coeffTable[i].y = wave(waveId).coeff.imag();

         //fill starId table
         starIdTable[i] = wave(waveId).starId;

         //invertflag table
         invertFlagTable[i] = stars_[starIdTable[i]].invertFlag;

      }

      IntVec<D> indiceBz;
      for(int i = 0; i < nStar_; i++) {
         for(int j = 0; j < D; j++) {
            indiceBz[j] = waveBz_[ i * D + j];
         }

         int rank = 0;
         //check if wave exists in dft
         if (!wave(indiceBz).implicit) {
            mesh_->shift(indiceBz);
            
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }

         }
         else {
            //wave does not exists in dft. have to find explicit pair
            mesh_->shift(indiceBz);
            for (int j = 0; j < D; ++j) {
               indiceBz[j] = -indiceBz[j];
            }
            mesh_->shift(indiceBz);

            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }

         }
         star2rank[i] = rank;
      }

      

      cudaMalloc((void**)&dksq, sizeof(cufftReal) * nStar_ * unitCellPtr_->nParams());
      //need to explicitly zero all vectors. 
      //I am surprised the code didnt explode the first time
      makeDksq(unitCell);
      

      //copy look up table to device
      //check cufftcomplex operator overload
      cudaMemcpy(rankTable_, rankTable, sizeof(int) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(partnerIdTable_, partnerIdTable, sizeof(int) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(implicitTable_, implicitTable, sizeof(bool) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(coeffTable_, coeffTable, sizeof(cufftComplex) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(starIdTable_, starIdTable, sizeof(int) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(invertFlagTable_, invertFlagTable, sizeof(int) * mesh_->size(), cudaMemcpyHostToDevice);
      cudaMemcpy(star2rank_, star2rank, sizeof(int) * mesh_->size(), cudaMemcpyHostToDevice);

      //delete host memory
      delete[] rankTable;
      delete[] partnerIdTable;
      delete[] implicitTable;
      delete[] coeffTable;
      delete[] starIdTable;
      delete[] invertFlagTable;
      delete[] star2rank;

   }

   template <int D>
   void Basis<D>::update(const UnitCell<D>& unitCell)
   {
      makeDksq(unitCell);

   }


#ifndef GPU_OUTER
   template <int D>
   void Basis<D>::convertFieldComponentsToDft(DArray<double>& components, RDFieldDft<D>& dft)
   {
      std::complex<double> coeff;
      double z;
      double z1;
      double z2;
      IntVec<D> G2;
      int partnerId;
      //loop through a factor of 2 more indices. 
      //Difficult to optimize without specialization or making new objects
      int kSize = 1;
      for(int i = 0; i < D; ++i) {
         if( i < D - 1) {
            kSize *= mesh_->dimension(i);
         } else {
            kSize *= (mesh_->dimension(i)/ 2 + 1);
         }
      }
      
      IntVec<D> offsets;
      offsets[D-1] = 1;
      for (int j = D-1; j > 0; --j) {
         if (j == D-1) {
            offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
         } else {
            offsets[j-1] = offsets[j]*mesh_->dimension(j);
         }
      }

      cufftComplex* temp = new cufftComplex[kSize];

      for (int i = 0; i < mesh_->size(); ++i) {
         
         //if last indice of waves is not > n3/2+1
         IntVec<D> waveId = mesh_->position(i);
         int starId = wave(waveId).starId;
         if (!wave(waveId).implicit) {
            coeff = wave(waveId).coeff;
            //determining the rank of RFieldDft
            
            int rank = 0;
            for (int j = 0; j < D; ++j) {
               rank += wave(waveId).indicesDft[j] * offsets[j];
            }
            //std::cout<<rank<<std::endl;

            switch (stars_[starId].invertFlag) {
               case 0 :
                  z = components[starId];
                  coeff = z*coeff;
                  //dft[rank][0] = coeff.real();
                  //dft[rank][1] = coeff.imag();
                  temp[rank].x = coeff.real();
                  temp[rank].y = coeff.imag();
                  break;
               case 1 :
                  z1 = components[starId];
                  for(int j = 0; j < D; ++j){
                     G2[j] = -waveId[j];
                  }
                  mesh_->shift(G2);
                  partnerId = mesh_->rank(G2);
                  z2 = components[partnerId];
                  /*//debug code
                  if(stars_[partnerId].invertFlag != -1) {
                     std::cout<<"This star is "<<starId<<" and my partner is "
                     <<partnerId<<std::endl;
                     std::cout<<"I have invertFlag of "<<stars_[starId].invertFlag
                     <<" and he has "<<stars_[partnerId].invertFlag<<std::endl;
                     std::cout<<" WaveBz "<<stars_[starId].waveBz <<" and he "
                     <<stars_[partnerId].waveBz<<std::endl;
                     std::cout<<"Wave dft indices"
                     <<waves_[stars_[starId].beginId].indicesDft<<"and he is"
                     <<waves_[stars_[partnerId].beginId].indicesDft<<std::endl;
                     std::cout<<"My previous wave has implicit of "
                     <<waves_[stars_[starId].beginId - 1].implicit<<std::endl;
                     std::cout<<"And position of "
                     <<waves_[stars_[starId].beginId-1].indicesDft<<std::endl;
                  }
                  //end debug*/
                  UTIL_CHECK(stars_[partnerId].invertFlag == -1);
                  coeff = std::complex<double>(z1,-z2)*coeff/sqrt(2);
                  //dft[rank][0] = coeff.real();
                  //dft[rank][1] = coeff.imag();
                  temp[rank].x = coeff.real();
                  temp[rank].y = coeff.imag();
                  break;
               case -1 :
                  z1 = components[starId];
                  for(int j = 0; j < D; ++j){
                     G2[j] = -waveId[j];
                  }
                  mesh_->shift(G2);
                  partnerId = mesh_->rank(G2);
                  z2 = components[partnerId];
                  UTIL_CHECK(stars_[partnerId].invertFlag == 1);
                  coeff = std::complex<double>(z2, z1) * coeff / sqrt(2);
                  //dft[rank][0] = coeff.real();
                  //dft[rank][1] = coeff.imag();
                  temp[rank].x = coeff.real();
                  temp[rank].y = coeff.imag();
                  break;
            }

         } else {
            //do nothing. Dft does not need allocation for this
         }
      }
      cudaMemcpy(dft.cDField(), temp, kSize * sizeof(cufftComplex), cudaMemcpyHostToDevice);
      delete[] temp;
   }
#else
      static __global__ void convertC2DftHelper(cufftReal* components, cufftComplex* dft, int* rankTable,
      int* partnerIdTable, bool* implicitTable, cufftComplex* coeffTable, int* starIdTable,
      int* invertFlagTable, int size) {

      int nThreads = blockDim.x * gridDim.x;
      int startID  = blockIdx.x * blockDim.x + threadIdx.x;

      cufftComplex coeff;
	  cufftComplex temp;
      double z;
      double z1;
      double z2;

      for(int i = startID; i < size; i += nThreads) {
         if(!implicitTable[i]) {
			 switch (invertFlagTable[i]) {
			 case 0:
				 z = components[starIdTable[i]];
				 coeff.x = z * coeffTable[i].x;
				 coeff.y = z * coeffTable[i].y;
				 dft[rankTable[i]] = coeff;
				 break;
			 case 1:
				 z1 = components[starIdTable[i]];
				 z2 = components[partnerIdTable[i]];
				 temp.x = z1;
				 temp.y = -z2;
				 //coeff = temp * coeffTable[i] / sqrt(2);
				 coeff = cuCmulf(temp, coeffTable[i]);
				 coeff.x /= sqrt(2.0);
				 coeff.y /= sqrt(2.0);
				 //cuFloatComplex cuCmulf(cuFloatComplex x,cuFloatComplex y)
				 //coeff.x = (temp.x * coeffTable[i].x - temp.y * coeffTable[i].y) / sqrt(2);
				 //coeff.y = (temp.x * coeffTable[i].y + temp.y * coeffTable[i].x) / sqrt(2);
				dft[rankTable[i]] = coeff;
				break;
            case -1:
               z1 = components[starIdTable[i]];
               z2 = components[partnerIdTable[i]];
               temp.x = z2;
               temp.y = z1;
			   //coeff = temp * coeffTable[i] / sqrt(2);
			   coeff = cuCmulf(temp, coeffTable[i]);
			   coeff.x /= sqrt(2.0);
			   coeff.y /= sqrt(2.0);
			   //coeff.x = (temp.x * coeffTable[i].x - temp.y * coeffTable[i].y) / sqrt(2);
			   //coeff.y = (temp.x * coeffTable[i].y + temp.y * coeffTable[i].x) / sqrt(2);
               dft[rankTable[i]] = coeff;
               break;
            }
         }
      }
   }
	template <int D>
	void Basis<D>::convertFieldComponentsToDft(RDField<D>& components, RDFieldDft<D>& dft)
	{
		
      //2 critic of this system. 
      // 1) pass a lot of table. Issue is not time but code management/readability
      // 2) some potential time save if pass function instead ( no need to loop again to make tables)
		convertC2DftHelper<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(components.cDField(), dft.cDField(), rankTable_,
                                          partnerIdTable_, implicitTable_, coeffTable_, 
                                          starIdTable_, invertFlagTable_, mesh_->size() );

	}
#endif

#ifndef GPU_OUTER
   template <int D>
   void Basis<D>::convertFieldDftToComponents(RDFieldDft<D>& dft, DArray<double>& components)
   {
      //if its not too distrupting maybe use consistent names for logical size?
      int nStars = nStar_;
      UTIL_CHECK(nStars == components.capacity());
      IntVec<D> indiceBz;
      std::complex<double> coeff;
      std::complex<double> temp;
      cufftComplex z;

      int kSize = 1;
      for(int i = 0; i < D; ++i) {
         if( i < D - 1) {
            kSize *= mesh_->dimension(i);
         } else {
            kSize *= (mesh_->dimension(i)/ 2 + 1);
         }
      }

      IntVec<D> offsets;
      offsets[D-1] = 1;
      for (int j = D-1; j > 0; --j) {
         if(j == D-1) {
            offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
         } else {
            offsets[j-1] = offsets[j]*mesh_->dimension(j);
         }
      }

      cufftComplex* kTemp = new cufftComplex[kSize];
      cudaMemcpy(kTemp, dft.cDField(), kSize * sizeof(cufftComplex), cudaMemcpyDeviceToHost);

      for (int i = 0; i < nStars; ++i) {
         for(int j = 0; j < D; ++j) {
            indiceBz[j] = waveBz_[i*D + j];
         }
         coeff = wave(indiceBz).coeff;

         //check if wave exists in dft
         if (!wave(indiceBz).implicit) {
            mesh_->shift(indiceBz);
            //std::cout<<"Yoohoo"<<std::endl;
            //might be good to rethink how elements of RFieldDft is accessed
            
            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            //std::cout<<offsets<<std::endl;
            //std::cout<<indiceBz<<std::endl;
            //std::cout<<rank<<std::endl;
            z.x = kTemp[rank].x;
            z.y = kTemp[rank].y;
            //Note: the code implies that we do not need indicesDft from Wave? 
         } else {
         //wave does not exists in dft. have to find explicit pair
            mesh_->shift(indiceBz);
            for (int j = 0; j < D; ++j) {
               indiceBz[j] = -indiceBz[j];
            }
            mesh_->shift(indiceBz);


            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1) {
                  offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
               } else {
                  offsets[j-1] = offsets[j]*mesh_->dimension(j);
               }
            }

            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            //std::cout<<"Yoohoo"<<std::endl;
            //std::cout<<rank<<std::endl;
            z.x = kTemp[rank].x;
            z.y = kTemp[rank].y;
            z.y = -z.y;
         }
         
         //std::cout<<i<<std::endl;
         //reintepret cast is not used since coding standards is old
         temp = std::complex<double>(z.x,z.y);
         //std::cout<<temp<<std::endl;
         //assign value to components
         switch(stars_[i].invertFlag) {
            case 0 :
               components[i] = (temp/coeff).real();
               break;
            case 1 :
               components[i] = (temp/coeff).real() * sqrt(2);
               break;
            case -1 :
               components[i] = (temp/coeff).imag() * sqrt(2);
               break;
         }
      }
      delete[] kTemp;
   }
#else

      static __global__ void convertDft2CHelper(cufftComplex* dft, cufftReal* components,
      bool* implicitTable, cufftComplex* coeffTable, int* invertFlagTable,
      int* star2rank, int size) {

      int nThreads = blockDim.x * gridDim.x;
      int startID  = blockIdx.x * blockDim.x + threadIdx.x;
      cufftComplex coeff;

      //figure out if coeffTable[i] is right
      cufftComplex z;
      for(int i = startID; i < size; i += nThreads) {

         coeff = coeffTable[i];
         if(!implicitTable[i]) {
            z.x = dft[star2rank[i]].x;
            z.y = dft[star2rank[i]].y;
         } else {
            z.x = dft[star2rank[i]].x;
            z.y = dft[star2rank[i]].y;
            z.y = -z.y;
         }
         switch(invertFlagTable[i]) {
         case 0:
            //components[i] = (z / coeff).x;
			components[i] = cuCdivf(z, coeff).x;
            break;
         case 1:
            //components[i] = (z / coeff).x * sqrt(2);
			components[i] = cuCdivf(z, coeff).x * sqrt(2.0);
            break;
         case -1:
            //compoenents[i] = (z / coeff).y * sqrt(2);
			 components[i] = cuCdivf(z, coeff).y * sqrt(2.0);
            break;
         }
      }
   }

	template <int D>
	void Basis<D>::convertFieldDftToComponents(RDFieldDft<D>& dft, RDField<D>& components)
	{
      //type something

		int nStars = nStar_;
		UTIL_CHECK(nStars == components.capacity());
		convertDft2CHelper<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(dft.cDField(), components.cDField(), implicitTable_,
         coeffTable_, invertFlagTable_, star2rank_, nStars);

      //cheating here
      //normally need to convert from starId to its corresponding waveBz
      //from waveBz convert to shift to be in mesh get the rank of that
      //and get its real position from waveId. Since the wave was not reorder
      // waveId and input indices is the same thing. waveId[itr.rank()] = itr.rank()
      // Do we now need the waveBz from starId and get its coressponding rank?

	}
#endif



   template <int D>
   int Basis<D>::nBasis() const
   {   
      int count = 0;
      for (int i = 0; i < stars_.capacity(); ++i) {
         if (!stars_[i].cancel) {
            count++;
         }   
      }   
      return count;
   }   

   template <int D>
   typename Basis<D>::Wave& Basis<D>::wave(IntVec<D> vector)
   {   
      if (!mesh_->isInMesh(vector)) {
         mesh_->shift(vector);
      }   
      return waves_[waveId_[mesh_->rank(vector)]];
   }   





   template <int D>
   void Basis<D>::makeDksq(const UnitCell<D>& unitCell)
   {
      
      unitCellPtr_ = &unitCell;
      
      /*
      // Initialize all elements to zero
      int i, j, p, q;
      //std::cout<<"unit cell ptr "<<unitCellPtr_->nParams()<<'\n';
      for (i = 0; i < unitCellPtr_->nParams(); ++i) {
         for (j=0; j < nStar_; ++j){
            dksq(i,j)=0.0;
         }
      }
      
      //let this be done on gpu
      //convert stars waveBz to a gpu array
      //find out time consumption
      int rank;
      for (i = 0; i < unitCellPtr_->nParams(); ++i) {
         for (j=0; j < nStar_; ++j){
            for (p=0; p < D; ++p){
               for (q=0; q < D; ++q){
                  rank = q + (p * D) + (i * (D * D));
                  dksq(i,j) = dksq(i,j) + (waveBz_[j*D + p] * waveBz_[j*D + q]
                                           * (double)unitCellPtr_->dkkBasis(rank) );
  
               }
            }
         }
      }
      */

      makeDksqHelper<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>
         (dksq, waveBz_d, unitCellPtr_->dkkBasis_d,
          unitCellPtr_->nParams(), nStar_, D);
   }


}
}

#endif
