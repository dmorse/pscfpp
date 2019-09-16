#ifndef PSSP_GPU_AM_ITERATOR_TPP
#define PSSP_GPU_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <pssp_gpu/System.h>
#include <util/format/Dbl.h>
//#include <Windows.h>
#include <thrust/count.h>
#include <thrust/reduce.h>
//#include <sys/time.h>

static __global__ void pointWiseAdd(cufftReal* result, const cufftReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] += rhs[i];
   }
}

static __global__ void subtractUniform(cufftReal* result, cufftReal rhs, int size) {
	int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = startID; i < size; i += nThreads) {
		result[i] -= rhs;
	}
}

static __global__ void pointWiseSubtract(cufftReal* result, const cufftReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] -= rhs[i];
   }
}

static __global__ void pointWiseBinarySubtract(const cufftReal* a, const cufftReal* b, cufftReal* result, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = a[i] - b[i];
   }
}

static __global__ void pointWiseBinaryAdd(const cufftReal* a, const cufftReal* b, cufftReal* result, int size) {
	int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = startID; i < size; i += nThreads) {
		result[i] = a[i] + b[i];
	}
}

static __global__ void pointWiseAddScale(cufftReal* result, const cufftReal* rhs,float scale, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] += scale * rhs[i];
   }
}

//the 1 is a placeholder for dr
static __global__ void AmIsConvergedHelper(cufftReal* out, int size) {
	int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;
	cufftReal temp;
	for (int i = startID; i < size; i += nThreads) {
		temp = (out[i] - 1) * (out[i] - 1) * 1;
		out[i] = temp;
	}
}

static __global__ void AmHelper(cufftReal* out, cufftReal* present, cufftReal* iPast, cufftReal* jPast, int size) {
	int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = startID; i < size; i += nThreads) {
		out[i] += (present[i] - iPast[i]) * (present[i] - jPast[i]);
	}
}

static __global__ void AmHelperVm(cufftReal* out, cufftReal* present, cufftReal* iPast, int size) {
	int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;
	for (int i = startID; i < size; i += nThreads) {
		out[i] += (present[i] - iPast[i]) * (present[i]);
	}
}

namespace Pscf {
namespace Pssp_gpu 
{

   using namespace Util;

   template <int D>
   AmIterator<D>::AmIterator()
    : Iterator<D>(),
      epsilon_(0),
      lambda_(0),
      nHist_(0),
      maxHist_(0)
   { setClassName("AmIterator"); }

   template <int D>
   AmIterator<D>::AmIterator(System<D>* system)
    : Iterator<D>(system),
      epsilon_(0),
      lambda_(0),
      nHist_(0),
      maxHist_(0)
   { setClassName("AmIterator"); }

   template <int D>
   AmIterator<D>::~AmIterator()
   {}

   template <int D>
   void AmIterator<D>::readParameters(std::istream& in)
   {
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
   }

   template <int D>
   void AmIterator<D>::allocate()
   {
      devHists_.allocate(maxHist_+1);
      omHists_.allocate(maxHist_+1);

      wArrays_.allocate(systemPtr_->mixture().nMonomer());
      dArrays_.allocate(systemPtr_->mixture().nMonomer());
      tempDev.allocate(systemPtr_->mixture().nMonomer());

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         wArrays_[i].allocate(systemPtr_->mesh().size());
         dArrays_[i].allocate(systemPtr_->mesh().size());
         tempDev[i].allocate(systemPtr_->mesh().size());
      }
	  
	  histMat_.allocate(maxHist_ + 1);


   }

   template <int D>
   int AmIterator<D>::solve()
   {
      //clock_t time_begin;
      //clock_t time_end;


      //solve the SCFT equations once
      //assumes basis.makeBasis() has been called
      //assumes AmIterator.allocate() has been called
      //consider making dft arrays size of one monomer and reuse the same array

      //check if setup?
      //systemPtr_->fft().setup(systemPtr_->wFieldGrid(0), systemPtr_->wFieldDft(0));
	  
      std::cout<<"Hahahahahahaha-----------";  
      
	   //requires you to read everything as wFieldGrids
      systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                    systemPtr_->cFieldGrids());

      systemPtr_->mixture().computeTStress(systemPtr_->basis());
                       for (int m=0; m<(systemPtr_->unitCell()).nParams() ; ++m){
                std::cout<<"Stress"<<m<<"\t"<<"="<< systemPtr_->mixture().TStress[m]<<"\n";
                 std::cout<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).params()[m]<<"\n";
              }  

	  //compute error at each mesh points
	  
      //check for convergence else resolve SCFT equations with new Fields
      for (int itr = 1; itr <= maxItr_; ++itr) {
         
         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr-1;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

		 SYSTEMTIME timeStart, timeEnd;
		 double timeElapsed;
		 GetSystemTime(&timeStart);

         computeDeviation();
         
		 GetSystemTime(&timeEnd);
		 timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
		 timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
		 timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
		 std::cout << " Time for computeDeviation ="
			 << Dbl(timeElapsed, 18, 11) << 's' << std::endl;

         std::cout<<"  Iteration  "<<itr<<std::endl;
		 
		 
		if(isConverged()) {
         //if (whatdo) {
            return 0;
         } else {
            //resize history based matrix appropriately
            //consider making these working space local

            /*gettimeofday(&tv, &tz);
            timeStart = (double)tv.tv_sec +
               (double)tv.tv_usec / 1000000.0;*/
			 GetSystemTime(&timeStart);

            if (itr <= maxHist_ + 1) {
               if (nHist_ > 0) {
                  invertMatrix_.allocate(nHist_, nHist_);
                  coeffs_.allocate(nHist_);
                  vM_.allocate(nHist_);
               }
            }
            
            minimizeCoeff(itr);
			            
			GetSystemTime(&timeEnd);
			timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			std::cout << " Time for minimizeCoeff ="
				<< Dbl(timeElapsed, 18, 11) << 's' << std::endl;
			
			//output old fields for debugging
			/*if (itr == 2) {
				std::ofstream outFile;
				std::string filename = "out/calc/omega2_old";
				systemPtr_->fileMaster().openOutputFile(filename, outFile);
				systemPtr_->writeRFields(outFile, systemPtr_->wFieldGrids());
				outFile.close();
			}*/

			GetSystemTime(&timeStart);
            buildOmega(itr);

			//cudaDeviceSynchronize();
			//takes about 1E-02 seconds if previous turned off
			GetSystemTime(&timeEnd);
			timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			std::cout << " Time for buildOmega ="
				<< Dbl(timeElapsed, 18, 11) << 's' << std::endl;
            

            /*gettimeofday(&tv, &tz);
            timeEnd = (double)tv.tv_sec + 
               (double)tv.tv_usec / 1000000.0;
			*/

            if (itr <= maxHist_) {
               //will deallocate when out of scope
               if (nHist_ > 0) {
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }

            
            /*for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               systemPtr_->basis().convertFieldComponentsToDft( 
                                    systemPtr_->wField(j),
                                    systemPtr_->wFieldDft(j));
               systemPtr_->fft().inverseTransform(systemPtr_->wFieldDft(j), 
                                                   systemPtr_->wFieldGrid(j));
            }
			*/

            /*gettimeofday(&tv, &tz);
            timeStart = (double)tv.tv_sec +
               (double)tv.tv_usec / 1000000.0;*/
			GetSystemTime(&timeStart);

            systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                          systemPtr_->cFieldGrids());
            /*gettimeofday(&tv, &tz);
            timeEnd = (double)tv.tv_sec + 
               (double)tv.tv_usec / 1000000.0;
            std::cout<<" Time for mixture compute ="
               << Dbl(timeEnd - timeStart,18,11)<<'s'<<std::endl;*/
			GetSystemTime(&timeEnd);
			timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			std::cout <<" Time for mixture compute ="
				<< Dbl(timeElapsed, 18, 11) << 's' << std::endl;

            
			/*GetSystemTime(&timeStart);
            for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
               systemPtr_->fft().forwardTransform(systemPtr_->cFieldGrid(i),
                                                  systemPtr_->cFieldDft(i));
               systemPtr_->basis().convertFieldDftToComponents(
                                    systemPtr_->cFieldDft(i),
                                    systemPtr_->cField(i));
            }
			GetSystemTime(&timeEnd);
			timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			std::cout << " Time for within loop transform ="
				<< Dbl(timeElapsed, 18, 11) << 's' << std::endl;

			*/
         }

      }

      //should not reach here. iterated more than maxItr. Not converged
      return 1;
   }

#ifndef GPU_OUTER
   template <int D>
   void AmIterator<D>::computeDeviation()
   {

      omHists_.append(systemPtr_->wFields());

      for (int i = 0 ; i < systemPtr_->mixture().nMonomer(); ++i) {
         for (int j = 0; j < systemPtr_->basis().nStar() - 1; ++j) {
            tempDev[i][j] = 0;
         }
      }
      //the form for this is slightly different for 3 species
      //Almost impossible to write good code here if using Interaction class
      //over ChiInteraction
      DArray<double> temp;
      temp.allocate(systemPtr_->basis().nStar() - 1);

      for(int i = 0; i < systemPtr_->basis().nStar() - 1; ++i) {
         temp[i] = 0;
      }
      for(int i = 0; i < systemPtr_->mixture().nMonomer() ++i) {
         for(int j = 0; j < systemPtr_->basis().nStar() -1; ++j) {
            temp[j] += systemPtr_->wField(i)[j + 1];
         }
      }
      
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
               tempDev[i][k] += systemPtr_->interaction().chi(i,j) *
                              systemPtr_->cField(j)[k + 1];
            }
         }

         for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
            tempDev[i][k] += ((temp[k] / systemPtr_->mixture().nMonomer())
                             - systemPtr_->wField(i)[k + 1]);
         }
      }

      //ultimately might be slow.copys the entire array. Better to just move
      //pointers
      devHists_.append(tempDev);

      //test code for IteratorTest.testComputeDeviation
      //should be all zero
      /*for(int i = 0; i < systemPtr_->mixture().nMonomer();i++){
         std::cout<<"THis is devfield of "<<i<<std::endl;
         for(int j = 0; j < systemPtr_->basis().nStar();j++){
            std::cout<<Dbl(devHists_[0][i][j])<<std::endl;
         }
      }*/
   }
#else
   template <int D>
   void AmIterator<D>::computeDeviation()
   {
	   
	   //average fields to zero
	   double average = 0;
	   for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		   thrust::device_ptr<cufftReal> avgPtr(systemPtr_->wFieldGrid(i).cDField());
		   average += thrust::reduce(avgPtr, avgPtr + systemPtr_->mesh().size());
	   }
	   average /= (systemPtr_->mixture().nMonomer() * systemPtr_->mesh().size());
	   for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		   subtractUniform << <1024, 256 >> > (systemPtr_->wFieldGrid(i).cDField(), average, systemPtr_->mesh().size());
	   }

      omHists_.append(systemPtr_->wFieldGrids());

      for (int i = 0 ; i < systemPtr_->mixture().nMonomer(); ++i){
         assignUniformReal<<<1024, 256>>>(tempDev[i].cDField(), 0, systemPtr_->mesh().size());
      }

      RDField<D> temp;
      temp.allocate(systemPtr_->mesh().size());
      assignUniformReal<<<1024, 256>>>(temp.cDField(), 0, systemPtr_->mesh().size());

	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  pointWiseAdd << <1024, 256 >> > (temp.cDField(), systemPtr_->wFieldGrid(i).cDField(), systemPtr_->mesh().size());
      }
      
	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
			  pointWiseAddScale << <1024, 256 >> > (tempDev[i].cDField(),
				  systemPtr_->cFieldGrid(j).cDField(),
				  systemPtr_->interaction().chi(i, j),
				  systemPtr_->mesh().size());
		  }

		  pointWiseAddScale << <1024, 256 >> > (tempDev[i].cDField(),
			  temp.cDField(), 1.0 / systemPtr_->mixture().nMonomer(),
			  systemPtr_->mesh().size());
	  }
	  //set new fields to average to zero
	  /*for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  thrust::device_ptr<cufftReal> avgPtr(tempDev[i].cDField());
		  average += thrust::reduce(avgPtr, avgPtr + systemPtr_->mesh().size());
		  std::cout << "this is predivided " << average << std::endl;
	  }
	  
	  average /= (systemPtr_->mixture().nMonomer() * systemPtr_->mesh().size());*/

	  //std::cout << "this is average: "<<average << std::endl;
	  //std::cout << " this is mesh size " << systemPtr_->mesh().size() << std::endl;
	  //std::cout << "this is monomer " << systemPtr_->mixture().nMonomer() << std::endl;
	  /*if (nHist_ == 0) {

		  std::ofstream outFile;
		  std::string filename = "out/calc/wbar";
		  systemPtr_->fileMaster().openOutputFile(filename, outFile);

		  cufftReal* tempArray;
		  tempArray = new cufftReal[systemPtr_->mesh().size()];
		  cudaMemcpy(tempArray, tempDev[0].cDField(),
			  systemPtr_->mesh().size() * sizeof(cufftReal), cudaMemcpyDeviceToHost);

		  for (int i = 0; i < systemPtr_->mesh().size(); i++) {
			  outFile << "  " << Dbl(tempArray[i], 18, 11) << std::endl;
		  }
		  outFile.close();
	  }*/
	  //again the average is known exactly here

	  average = 0;
	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  thrust::device_ptr<cufftReal> avgPtr(tempDev[i].cDField());
		  average += thrust::reduce(avgPtr, avgPtr + systemPtr_->mesh().size());
	  }
	  average /= (systemPtr_->mixture().nMonomer() * systemPtr_->mesh().size());
	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  subtractUniform << <1024, 256 >> > (tempDev[i].cDField(), average, systemPtr_->mesh().size());
	  }

	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
		  pointWiseSubtract << <1024, 256 >> >(tempDev[i].cDField(),
			  systemPtr_->wFieldGrid(i).cDField(),
			  systemPtr_->mesh().size());
	  }
         
      

      //ultimately might be slow.copys the entire array. Better to just move
      //pointers
      devHists_.append(tempDev);
   }
#endif
   
#ifndef GPU_OUTER
   template <int D>
   bool AmIterator<D>::isConverged()
   {
      double error;
      double dError = 0;
      double wError = 0;

      for ( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         for ( int j = 0; j < systemPtr_->basis().nStar() - 1; j++) {
            dError += devHists_[0][i][j] * devHists_[0][i][j];

            //the extra shift is due to the zero indice coefficient being
            //exactly known
            wError += systemPtr_->wField(i)[j+1] * systemPtr_->wField(i)[j+1];
         }
      }
      std::cout<<" dError------ :"<<Dbl(dError)<<std::endl;
      std::cout<<" wError :"<<Dbl(wError)<<std::endl;
      error = sqrt(dError / wError);
      std::cout<<"  Error  :"<<Dbl(error)<<std::endl;
      if (error < epsilon_) {
         return true;
      } else {
         return false;
      }
   }
#else
   template <int D>
   bool AmIterator<D>::isConverged()
   {
#if 0
      double error;

	  assignReal << <1024, 256 >> > (tempDev[0].cDField(), systemPtr_->cFieldGrid(0).cDField(), systemPtr_->mesh().size());
	  for (int i = 1; i < systemPtr_->mixture().nMonomer(); i++) {
		  pointWiseAdd << <1024, 256 >> > (tempDev[0].cDField(), systemPtr_->cFieldGrid(i).cDField(), systemPtr_->mesh().size());
	  }
	  //take phi(r) and make (phi(r) - 1)^2 * dr
	  AmIsConvergedHelper << <1024, 256 >> > (tempDev[0].cDField(), systemPtr_->mesh().size());

	  //sum (phi(r) - 1)^2 dr
	  thrust::device_ptr<cufftReal> errorPtr(tempDev[0].cDField());
	  error = thrust::reduce(errorPtr, errorPtr + systemPtr_->mesh().size());

	  //1/v * sum(phi(r) - 1)^2 dr
	  error /= 1000000;

      std::cout<<"  Error  :"<<Dbl(error)<<std::endl;
      if (error < epsilon_) {
         return true;
      } else {
         return false;
      }
#endif
	  
	  double error;
	  cufftReal dError = 0;
	  cufftReal wError = 0;

	  for (int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
		  //squaring the devHists, store in temp work array
		  pointwiseMul << <1024, 256 >> >(devHists_[0][i].cDField(),
			  devHists_[0][i].cDField(), tempDev[i].cDField(), systemPtr_->mesh().size());
		  thrust::device_ptr<cufftReal> dErrorPtr(tempDev[i].cDField());
		  //reducing the squared value
		  dError += thrust::reduce(dErrorPtr, dErrorPtr + systemPtr_->mesh().size());

		  //squaring wError
		  pointwiseMul << <1024, 256 >> >(systemPtr_->wFieldGrid(i).cDField(),
			  systemPtr_->wFieldGrid(i).cDField(), tempDev[i].cDField(),
			  systemPtr_->mesh().size());
		  thrust::device_ptr<cufftReal> wErrorPtr(tempDev[i].cDField());
		  wError += thrust::reduce(wErrorPtr, wErrorPtr + systemPtr_->mesh().size());

	  }
	  std::cout << " dError ---------:" << Dbl(dError) << std::endl;
	  std::cout << " wError :" << Dbl(wError) << std::endl;
	  error = sqrt(dError / wError);
	  std::cout << "  Error  :" << Dbl(error) << std::endl;
	  if (error < epsilon_) {
		  return true;
	  }
	  else {
		  return false;
	  }
   }
#endif


#ifndef GPU_OUTER
   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
      //clock_t time_begin;
      //clock_t time_end;
      if (itr == 1) {
         //do nothing
      } else {

         double elm;

         for (int i = 0; i < nHist_; ++i) {
            for (int j = i; j < nHist_; ++j) {


               invertMatrix_(i,j) = 0;

               for (int k = 0; k < systemPtr_->mixture().nMonomer(); ++k) {
                  elm = 0;
                  for (int l = 0; l < systemPtr_->basis().nStar() - 1; ++l) {
                     elm +=
                        (  (devHists_[0][k][l] - devHists_[i+1][k][l]) *
                           (devHists_[0][k][l] - devHists_[j+1][k][l]) );
                  }

                  invertMatrix_(i,j) += elm;


               }

               invertMatrix_(j,i) = invertMatrix_(i,j);
            }

            vM_[i] = 0;
            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
                  vM_[i] += ( (devHists_[0][j][k] - devHists_[i+1][j][k]) *
                               devHists_[0][j][k] );
               }
            }
         }

         if (itr == 2) {
            coeffs_[0] = vM_[0] / invertMatrix_(0,0);
         } else {

			 /*SYSTEMTIME timeStart, timeEnd;
			 double timeElapsed;
			 GetSystemTime(&timeStart);
			 */
            //time_begin = clock();
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);
            solver.solve(vM_, coeffs_);

			/*GetSystemTime(&timeEnd);
			timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			std::cout << " Time for inverting the matrix ="
				<< Dbl(timeElapsed, 18, 11) << 's' << std::endl;
			*/
            //time_end = clock();
            //std::cout<<" nHist_ is "<<nHist_<<std::endl;
            //std::cout<<" Time for LUSolver ="
            //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;
         }       
      }
   }
#else

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
	   SYSTEMTIME timeStart, timeEnd;
	   double timeElapsed;
      //clock_t time_begin;
      //clock_t time_end;
      if (itr == 1) {
         //do nothing
      } else {

         float elm;

		 GetSystemTime(&timeStart);
		 //clear last column and shift everything downwards if necessary
		 histMat_.clearColumn(nHist_);
		 //calculate the new values for d(k)d(k') matrix
		 for (int i = 0; i < maxHist_ + 1; ++i) {
			 if (i < nHist_ + 1) {
				 elm = 0;
				 for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
					 pointwiseMul << <1024, 256 >> > (devHists_[0][j].cDField(), 
						 devHists_[i][j].cDField(), tempDev[0].cDField(), 
						 systemPtr_->mesh().size());
					 thrust::device_ptr<cufftReal> elmPtr(tempDev[0].cDField());
					 elm += thrust::reduce(elmPtr, elmPtr + systemPtr_->mesh().size());
				 }
				 histMat_.evaluate(elm, nHist_, i);
			 }
		 }

		 //build new Umn array
		 for (int i = 0; i < nHist_; ++i) {
			 for (int j = i; j < nHist_; ++j) {
				 invertMatrix_(i, j) = histMat_.makeUmn(i, j, nHist_);
				 invertMatrix_(j, i) = invertMatrix_(i, j);
			 }
		 }

		 
		 GetSystemTime(&timeEnd);
		 timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
		 timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
		 timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
		 std::cout << " Time for forming elm ="
			 << Dbl(timeElapsed, 18, 11) << 's' << std::endl;

		 GetSystemTime(&timeStart);
		 for (int i = 0; i < nHist_; ++i) {
			 vM_[i] = histMat_.makeVm(i, nHist_);
		 }
		 
		 GetSystemTime(&timeEnd);
		 timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
		 timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
		 timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
		 std::cout << " Time for forming Vm ="
			 << Dbl(timeElapsed, 18, 11) << 's' << std::endl;


         if (itr == 2) {
            coeffs_[0] = vM_[0] / invertMatrix_(0,0);
         } else {

          
          GetSystemTime(&timeStart);
          
            
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);
            solver.solve(vM_, coeffs_);

         GetSystemTime(&timeEnd);
         timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
         timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
         timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
         std::cout << " Time for inverting the matrix ="
            << Dbl(timeElapsed, 18, 11) << 's' << std::endl;
         
			//debugging code
			 /*if (itr == 10) {
				 for (int i = 0; i < nHist_; ++i) {
					 std::cout << coeffs_[i] << std::endl;
				 }
			 }*/
         }       
      }
   }
#endif

#if 0

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
	   SYSTEMTIME timeStart, timeEnd;
	   double timeElapsed;
	   //clock_t time_begin;
	   //clock_t time_end;
	   if (itr == 1) {
		   //do nothing
	   }
	   else {

		   double elm;

		   GetSystemTime(&timeStart);
		   for (int i = 0; i < nHist_; ++i) {
			   for (int j = i; j < nHist_; ++j) {


				   invertMatrix_(i, j) = 0;

				   assignUniformReal << <1024, 256 >> > (tempDev[0].cDField(), 0.0, systemPtr_->basis().nStar() - 1);
				   for (int k = 0; k < systemPtr_->mixture().nMonomer(); ++k) {
					   AmHelper << <1024, 256 >> > (tempDev[0].cDField(), devHists_[0][k].cDField(),
						   devHists_[i + 1][k].cDField(), devHists_[j + 1][k].cDField(), systemPtr_->basis().nStar() - 1);
					   /*pointWiseBinarySubtract << <1024, 256 >> > (devHists_[0][k].cDField(),
					   devHists_[i + 1][k].cDField(), tempDev[0].cDField(),
					   systemPtr_->basis().nStar() - 1);
					   pointWiseBinarySubtract << <1024, 256 >> > (devHists_[0][k].cDField(),
					   devHists_[j + 1][k].cDField(), tempDev[1].cDField(),
					   systemPtr_->basis().nStar() - 1);
					   inPlacePointwiseMul << <1024, 256 >> > (tempDev[0].cDField(), tempDev[1].cDField(),
					   systemPtr_->basis().nStar() - 1);*/
				   }
				   thrust::device_ptr<cufftReal> elmPtr(tempDev[0].cDField());
				   invertMatrix_(i, j) = thrust::reduce(elmPtr, elmPtr + systemPtr_->basis().nStar() - 1);

				   invertMatrix_(j, i) = invertMatrix_(i, j);
			   }
		   }
		   GetSystemTime(&timeEnd);
		   timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
		   timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
		   timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
		   std::cout << " Time for forming elm ="
			   << Dbl(timeElapsed, 18, 11) << 's' << std::endl;

		   GetSystemTime(&timeStart);
		   for (int i = 0; i < nHist_; ++i) {
			   vM_[i] = 0;

			   assignUniformReal << <1024, 256 >> > (tempDev[0].cDField(), 0.0, systemPtr_->basis().nStar() - 1);
			   for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
				   AmHelperVm << <1024, 256 >> > (tempDev[0].cDField(), devHists_[0][j].cDField(),
					   devHists_[i + 1][j].cDField(), systemPtr_->basis().nStar() - 1);
				   /*pointWiseBinarySubtract<<<1024, 256>>>(devHists_[0][j].cDField(),
				   devHists_[i+1][j].cDField(), tempDev[0].cDField(),
				   systemPtr_->basis().nStar()-1);
				   inPlacePointwiseMul<<<1024, 256>>>(tempDev[0].cDField(), devHists_[0][j].cDField(),
				   systemPtr_->basis().nStar() - 1);*/
				   //reducing the squared value
			   }
			   //cudaDeviceSynchronize();
			   thrust::device_ptr<cufftReal> vmPtr(tempDev[0].cDField());
			   vM_[i] = thrust::reduce(vmPtr, vmPtr + systemPtr_->basis().nStar() - 1);
		   }
		   GetSystemTime(&timeEnd);
		   timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
		   timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
		   timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
		   std::cout << " Time for forming Vm ="
			   << Dbl(timeElapsed, 18, 11) << 's' << std::endl;


		   if (itr == 2) {
			   coeffs_[0] = vM_[0] / invertMatrix_(0, 0);
		   }
		   else {


			   GetSystemTime(&timeStart);

			   //time_begin = clock();
			   LuSolver solver;
			   solver.allocate(nHist_);
			   solver.computeLU(invertMatrix_);
			   solver.solve(vM_, coeffs_);

			   GetSystemTime(&timeEnd);
			   timeElapsed = (timeEnd.wMinute - timeStart.wMinute) * 60;
			   timeElapsed += (timeEnd.wSecond - timeStart.wSecond);
			   timeElapsed += (timeEnd.wMilliseconds - timeStart.wMilliseconds) / 1000.0;
			   std::cout << " Time for inverting the matrix ="
				   << Dbl(timeElapsed, 18, 11) << 's' << std::endl;

			   //time_end = clock();
			   //std::cout<<" nHist_ is "<<nHist_<<std::endl;
			   //std::cout<<" Time for LUSolver ="
			   //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;
		   }
	   }
   }
#endif

#ifndef GPU_OUTER
   template <int D>
   void AmIterator<D>::buildOmega(int itr)
   {

      if (itr == 1) {
         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            for (int j = 0; j < systemPtr_->basis().nStar() - 1; ++j) {
               systemPtr_->wField(i)[j+1] = omHists_[0][i][j+1] + 
                                             lambda_*devHists_[0][i][j];
            }
         }
      } else {
         //should be strictly correct. coeffs_ is a vector of size 1 if itr ==2

         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
               //extra shift in wArrays because omHists stores the first star
               wArrays_[j][k] = omHists_[0][j][k + 1];
               dArrays_[j][k] = devHists_[0][j][k];
            }
         }

         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
                  wArrays_[j][k] += coeffs_[i] * ( omHists_[i+1][j][k+1] - 
                                                   omHists_[0][j][k+1] );
                  dArrays_[j][k] += coeffs_[i] * ( devHists_[i+1][j][k] - 
                                                   devHists_[0][j][k] );
               }
            }
         }


         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            for (int j = 0; j < systemPtr_->basis().nStar() - 1; ++j) {
              systemPtr_->wField(i)[j+1] = wArrays_[i][j] + lambda_ * dArrays_[i][j];
               //std::cout<<systemPtr_->wField(i)[j+1]<<std::endl;

            }
         }

      }
   }
#else
   template <int D>
   void AmIterator<D>::buildOmega(int itr)
   {

      if (itr == 1) {
         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            assignReal<<<1024, 256>>>(systemPtr_->wFieldGrid(i).cDField(),
               omHists_[0][i].cDField(), systemPtr_->mesh().size());
            pointWiseAddScale<<<1024, 256>>>(systemPtr_->wFieldGrid(i).cDField(),
               devHists_[0][i].cDField(), lambda_ , systemPtr_->mesh().size());
         }
      } else {
         //should be strictly correct. coeffs_ is a vector of size 1 if itr ==2

         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            assignReal<<<1024, 256>>>(wArrays_[j].cDField(),
               omHists_[0][j].cDField(), systemPtr_->mesh().size());
            assignReal<<<1024, 256>>>(dArrays_[j].cDField(),
               devHists_[0][j].cDField(), systemPtr_->mesh().size());

         }

         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               //wArrays
               pointWiseBinarySubtract<<<1024, 256>>>(omHists_[i+1][j].cDField(),
                  omHists_[0][j].cDField(), tempDev[0].cDField(),
                   systemPtr_->mesh().size());
               pointWiseAddScale<<<1024, 256>>>(wArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i] , systemPtr_->mesh().size());

               //dArrays
               pointWiseBinarySubtract<<<1024, 256>>>(devHists_[i+1][j].cDField(),
                  devHists_[0][j].cDField(), tempDev[0].cDField(),
                  systemPtr_->mesh().size());
               pointWiseAddScale<<<1024, 256>>>(dArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i] , systemPtr_->mesh().size());
            }
         }


         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            assignReal<<<1024, 256>>>(systemPtr_->wFieldGrid(i).cDField(),
               wArrays_[i].cDField(), systemPtr_->mesh().size());
            pointWiseAddScale<<<1024, 256>>>(systemPtr_->wFieldGrid(i).cDField(),
               dArrays_[i].cDField(), lambda_ , systemPtr_->mesh().size());
         }

      }
   }
#endif
}
}

#endif
