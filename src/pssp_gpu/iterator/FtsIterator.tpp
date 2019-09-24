#ifndef PSSP_GPU_FTS_ITERATOR_TPP
#define PSSP_GPU_FTS_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FtsIterator.h"
#include <pssp_gpu/System.h>
#include <util/format/Dbl.h>
//#include <Windows.h>
#include <pssp_gpu/GpuResources.h>
#include <util/containers/FArray.h>
#include <sys/time.h>


namespace Pscf {
namespace Pssp_gpu {

   using namespace Util;

   template <int D>
   FtsIterator<D>::FtsIterator()
      : Iterator<D>(),
        epsilon_(0),
        lambda_(0),
        nHist_(0),
        maxHist_(0),
        cell_(0)
   {
      //temporary for testing
      setClassName("AmIterator");
   }

   template <int D>
   FtsIterator<D>::FtsIterator(System<D>* system)
      : Iterator<D>(system),
        epsilon_(0),
        lambda_(0),
        nHist_(0),
        maxHist_(0),
        cell_(0)
   {
      setClassName("AmIterator");
   }

   template <int D>
   FtsIterator<D>::~FtsIterator()
   {
      delete[] temp_;
      cudaFree(d_temp_);
   }

   template <int D>
   void FtsIterator<D>::readParameters(std::istream& in)
   {   
      cell_ = 0; //default value (fixed cell)
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "domain", cell_); 
   }

   template <int D>
   void FtsIterator<D>::allocate()
   {
      devHists_.allocate(maxHist_ + 1);
      omHists_.allocate(maxHist_ + 1);

      if (cell_){
         devCpHists_.allocate(maxHist_+1);
         CpHists_.allocate(maxHist_+1);
      }

      wArrays_.allocate(systemPtr_->mixture().nMonomer());
      dArrays_.allocate(systemPtr_->mixture().nMonomer());
      tempDev.allocate(systemPtr_->mixture().nMonomer());

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
          wArrays_[i].allocate(systemPtr_->mesh().size());
          dArrays_[i].allocate(systemPtr_->mesh().size());
          tempDev[i].allocate(systemPtr_->mesh().size());
      }
      
      histMat_.allocate(maxHist_ + 1);
      //allocate d_temp_ here i suppose
      cudaMalloc((void**)&d_temp_, NUMBER_OF_BLOCKS * sizeof(cufftReal));
      temp_ = new cufftReal[NUMBER_OF_BLOCKS];
   }

   template <int D>
   int FtsIterator<D>::solve()
   {
      
      struct timeval timeStart, timeEnd;
      struct timezone tz;

      gettimeofday(&timeStart, &tz);
      //data read in r-space
      systemPtr_->mixture().compute(systemPtr_->wFieldGrids(),
         systemPtr_->cFieldGrids());

      gettimeofday(&timeEnd, &tz);
      Log::file() <<"MDE time : "
                  <<Dbl( ( (double)(timeEnd.tv_sec - timeStart.tv_sec) + 
                           (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0 )
                         , 18, 11) <<'\n';
      //what if cell_ is not read?
      //why do you point to basis functions?
      gettimeofday(&timeStart, &tz);
      systemPtr_->mixture().computeTStress(systemPtr_->basis());
      gettimeofday(&timeEnd, &tz);
      Log::file() <<"Stress time : "
                  <<Dbl( ( (double)(timeEnd.tv_sec - timeStart.tv_sec) + 
                           (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0 )
                         , 18, 11) <<'\n';
     
      //for loop formatting all wrong
      //why is there an extra parenthesis in the pointer
      for (int m = 0; m< systemPtr_->unitCell().nParams() ; ++m){
         std::cout<<"Stress"<<m<<"\t"<<"="<< systemPtr_->mixture().TStress[m]<<"\n";
         std::cout<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).params()[m]<<"\n";
      } 


      //compute error at each mesh points

      //check for convergence else resolve SCFT equations with new Fields
      int itr;
      for (itr = 1; itr <= maxItr_; ++itr) {
        
         std::cout<< "iterations : " << itr << std::endl;
         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr - 1;
         }
         else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         computeDeviation();

         if (isConverged()) {
            if (itr > maxHist_ + 1) {
               invertMatrix_.deallocate();
               coeffs_.deallocate();
               vM_.deallocate();
            }
            //the space are misaligned
            //is it rational to print outputs here
            std::cout<<"----------CONVERGED----------"<< std::endl;
            if(cell_) {
               for (int m = 0; m < systemPtr_->unitCell().nParams() ; ++m){
                  std::cout<<"Stress"<<m<<"\t"<<"="<< systemPtr_->mixture().TStress[m]<<"\n";
                  std::cout<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).params()[m]<<"\n";
               }
            }
            return 0;
         }
         else {
            //resize history based matrix appropriately
            //consider making these working space local

            //GetSystemTime(&timeStart);


             if (itr <= maxHist_ + 1) {
                if (nHist_ > 0) {
                   invertMatrix_.allocate(nHist_, nHist_);
                   coeffs_.allocate(nHist_);
                   vM_.allocate(nHist_);
                }
             }


            int status = minimizeCoeff(itr);

            if (status == 1) {
               //abort the calculations and treat as failure (time out)
               //perform some clean up stuff
               invertMatrix_.deallocate();
               coeffs_.deallocate();
               vM_.deallocate();
               return 1;

            }


            buildOmega(itr);

            if (itr <= maxHist_) {
               if (nHist_ > 0) {
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }

            //GetSystemTime(&timeStart);
            gettimeofday(&timeStart, &tz);     
            systemPtr_->mixture().compute(systemPtr_->wFieldGrids(),
               systemPtr_->cFieldGrids());
         
            gettimeofday(&timeEnd, &tz);
            Log::file() <<"MDE time : "
                        <<Dbl( ( (double)(timeEnd.tv_sec - timeStart.tv_sec) + 
                                 (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0 )
                               , 18, 11) <<'\n';
     
            gettimeofday(&timeStart, &tz);      
            if (cell_){
               systemPtr_->mixture().computeTStress(systemPtr_->basis());
               
               for (int m = 0; m < systemPtr_->unitCell().nParams() ; ++m){
                  std::cout<<"Stress"<<m<<"\t"<<"="
                           << std::setprecision (15)
                           << systemPtr_->mixture().TStress[m]<<"\n";
               }
            }
            gettimeofday(&timeEnd, &tz);
            Log::file() <<"Stress time : "
                        <<Dbl( ( (double)(timeEnd.tv_sec - timeStart.tv_sec) + 
                                 (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0 )
                               , 18, 11) <<'\n';
            

   

         }
         
      }
      if (itr > maxHist_ + 1) {
         invertMatrix_.deallocate();
         coeffs_.deallocate();
         vM_.deallocate();
      }
      
      //should not reach here. iterated more than maxItr. Not converged
      return 1;
   }

   template <int D>
   void FtsIterator<D>::computeDeviation()
   {

      //need to average
      float average = 0;
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         average += reductionH(systemPtr_->wFieldGrid(i), systemPtr_->mesh().size());
      }
      average /= (systemPtr_->mixture().nMonomer() * systemPtr_->mesh().size());
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         subtractUniform << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (systemPtr_->wFieldGrid(i).cDField(), average, systemPtr_->mesh().size());
      }

      omHists_.append(systemPtr_->wFieldGrids());

      if (cell_) {
         CpHists_.append(systemPtr_->unitCell().params());
      }

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         assignUniformReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(tempDev[i].cDField(), 0, systemPtr_->mesh().size());
      }
      
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (tempDev[i].cDField(),
                                                                            systemPtr_->cFieldGrid(j).cDField(),
                                                                            systemPtr_->interaction().chi(i, j),
                                                                            systemPtr_->mesh().size());
            //this is a good add but i dont necessarily know if its right
            pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> > (tempDev[i].cDField(),
                                                                            systemPtr_->wFieldGrid(j).cDField(),
                                                                            -systemPtr_->interaction().indemp(i, j),
                                                                            systemPtr_->mesh().size()); 

         }
         
      }

      float sum_chi_inv = (float) systemPtr_->interaction().sum_inv;

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         
         pointWiseSubtractFloat << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(tempDev[i].cDField(),
                                                                             1/sum_chi_inv, 
                                                                             systemPtr_->mesh().size());
      }



      devHists_.append(tempDev);


      if (cell_){
         FArray<double, 6> tempCp;
         for (int i = 0; i<(systemPtr_->unitCell()).nParams(); i++){
            //format????
            tempCp [i] = -((systemPtr_->mixture()).TStress [i]);
         }
         devCpHists_.append(tempCp);
      }
   }

   template <int D>
   bool FtsIterator<D>::isConverged()
   {
      double error;
      double dError = 0;
      double wError = 0;
      //float temp = 0;
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         dError += innerProduct(devHists_[0][i], devHists_[0][i], systemPtr_->mesh().size());
         wError += innerProduct(systemPtr_->wFieldGrid(i), systemPtr_->wFieldGrid(i), systemPtr_->mesh().size());
      }

      if (cell_){
         for ( int i = 0; i < systemPtr_->unitCell().nParams(); i++) {
            dError +=  devCpHists_[0][i] *  devCpHists_[0][i];
            wError +=  systemPtr_->unitCell().params() [i] * systemPtr_->unitCell().params() [i];
         }
      }

      std::cout << " dError :" << Dbl(dError) << std::endl;
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

   template <int D>
   int FtsIterator<D>::minimizeCoeff(int itr)
   {
      if (itr == 1) {
         //do nothing
         histMat_.reset();
         return 0;
      }
      else {


         float elm, elm_cp;
         //clear last column and shift everything downwards if necessary
         histMat_.clearColumn(nHist_);
         //calculate the new values for d(k)d(k') matrix
         for (int i = 0; i < maxHist_ + 1; ++i) {
            if (i < nHist_ + 1) {

               elm = 0;
               for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
                  elm += FtsIterator<D>::innerProduct(devHists_[0][j], devHists_[i][j], systemPtr_->mesh().size());
               }
               histMat_.evaluate(elm, nHist_, i);
            }
         }

         //build new Umn array
         for (int i = 0; i < nHist_; ++i) {
            for (int j = i; j < nHist_; ++j) {
               invertMatrix_(i, j) = histMat_.makeUmn(i, j, nHist_);

               if (cell_){
                  elm_cp = 0;
                  for (int m = 0; m < systemPtr_->unitCell().nParams(); ++m){
                     elm_cp += ((devCpHists_[0][m] - devCpHists_[i+1][m]) * 
                                (devCpHists_[0][m] - devCpHists_[j+1][m])); 
                  }
                  invertMatrix_(i,j) += elm_cp;
               } 
               
               invertMatrix_(j, i) = invertMatrix_(i, j);
            }
         }

         for (int i = 0; i < nHist_; ++i) {
            vM_[i] = histMat_.makeVm(i, nHist_);

            if (cell_){
               elm_cp = 0;
               for (int m = 0; m < systemPtr_->unitCell().nParams(); ++m){
                  vM_[i] += ((devCpHists_[0][m] - devCpHists_[i+1][m]) *
                             (devCpHists_[0][m]));
               }
            }
         }

         if (itr == 2) {
            coeffs_[0] = vM_[0] / invertMatrix_(0, 0);
         }
         else {

            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);

            /*
            int status = solver.solve(vM_, coeffs_);
            if (status) {
               if (status == 1) {
                  //matrix is singular do something
                  return 1;
               }
               }*/
            solver.solve(vM_, coeffs_);
            //for the sake of simplicity during porting
            //we leaves out checks for singular matrix here
            //--GK 09 11 2019

         }
         return 0;
      }
   }

   template <int D>
   cufftReal FtsIterator<D>::innerProduct(const RDField<D>& a, const RDField<D>& b, int size) {

     switch(THREADS_PER_BLOCK){
     case 512:
       deviceInnerProduct<512><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 256:
       deviceInnerProduct<256><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 128:
       deviceInnerProduct<128><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 64:
       deviceInnerProduct<64><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 32:
       deviceInnerProduct<32><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 16:
       deviceInnerProduct<16><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 8:
       deviceInnerProduct<8><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 4:
       deviceInnerProduct<4><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 2:
       deviceInnerProduct<2><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     case 1:
       deviceInnerProduct<1><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cufftReal)>>>(d_temp_, a.cDField(), b.cDField(), size);
       break;
     }
      cudaMemcpy(temp_, d_temp_, NUMBER_OF_BLOCKS * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      cufftReal final = 0;
      cufftReal c = 0;
      //use kahan summation to reduce error
      for (int i = 0; i < NUMBER_OF_BLOCKS; ++i) {
         cufftReal y = temp_[i] - c;
         cufftReal t = final + y;
         c = (t - final) - y;
         final = t;
         
      }
      
      return final;
   }


   template<int D>
   cufftReal FtsIterator<D>::reductionH(const RDField<D>& a, int size) {
     reduction <<< NUMBER_OF_BLOCKS/2 , THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cufftReal) >> > (d_temp_, a.cDField(), size);
      cudaMemcpy(temp_, d_temp_, NUMBER_OF_BLOCKS/2  * sizeof(cufftReal), cudaMemcpyDeviceToHost);
      cufftReal final = 0;
      cufftReal c = 0;
      for (int i = 0; i < NUMBER_OF_BLOCKS/2 ; ++i) {
         cufftReal y = temp_[i] - c;
         cufftReal t = final + y;
         c = (t - final) - y;
         final = t;
      }
      return final;
   }

   template <int D>
   void FtsIterator<D>::buildOmega(int itr)
   {

      if (itr == 1) {
         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            assignReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(systemPtr_->wFieldGrid(i).cDField(),
            omHists_[0][i].cDField(), systemPtr_->mesh().size());
            pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(systemPtr_->wFieldGrid(i).cDField(),
            devHists_[0][i].cDField(), lambda_, systemPtr_->mesh().size());
         }

         if (cell_){
            for (int m = 0; m < (systemPtr_->unitCell()).nParams() ; ++m){
               systemPtr_->unitCell().SetParams( CpHists_[0][m] +lambda_* devCpHists_[0][m] ,m);
            }
            
            (systemPtr_->unitCell()).setLattice();
            systemPtr_->mixture().setupUnitCell(systemPtr_->unitCell());
            systemPtr_->basis().update(systemPtr_->unitCell());
            
            for (int m = 0; m < systemPtr_->unitCell().nParams(); ++m){
               std::cout<<"Parameter"<<m<<"\t"<<"="<<systemPtr_->unitCell().params()[m]<<"\n";
            }
         }

      }
      else {
         //should be strictly correct. coeffs_ is a vector of size 1 if itr ==2

         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            assignReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(wArrays_[j].cDField(),
               omHists_[0][j].cDField(), systemPtr_->mesh().size());
            assignReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(dArrays_[j].cDField(),
               devHists_[0][j].cDField(), systemPtr_->mesh().size());

         }

         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               //wArrays
               pointWiseBinarySubtract << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(omHists_[i + 1][j].cDField(),
                  omHists_[0][j].cDField(), tempDev[0].cDField(),
                  systemPtr_->mesh().size());
               pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(wArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i], systemPtr_->mesh().size());

               //dArrays
               pointWiseBinarySubtract << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(devHists_[i + 1][j].cDField(),
                  devHists_[0][j].cDField(), tempDev[0].cDField(),
                  systemPtr_->mesh().size());
               pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(dArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i], systemPtr_->mesh().size());
            }
         }
         
         for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
            assignReal << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(systemPtr_->wFieldGrid(i).cDField(),
               wArrays_[i].cDField(), systemPtr_->mesh().size());
            pointWiseAddScale << <NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >> >(systemPtr_->wFieldGrid(i).cDField(),
               dArrays_[i].cDField(), lambda_, systemPtr_->mesh().size());
         }



         if(cell_){
            
            for (int m = 0; m < systemPtr_->unitCell().nParams() ; ++m){
               wCpArrays_[m] = CpHists_[0][m];
               dCpArrays_[m] = devCpHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < systemPtr_->unitCell().nParams() ; ++m) {
                  wCpArrays_[m] += coeffs_[i] * ( CpHists_[i+1][m]-
                                                  CpHists_[0][m]);
                  dCpArrays_[m] += coeffs_[i] * ( devCpHists_[i+1][m]-
                                                  devCpHists_[0][m]);
               }
            } 

            struct timeval timeStart, timeEnd;
            struct timezone tz;

            gettimeofday(&timeStart, &tz);                
            for (int m = 0; m < systemPtr_->unitCell().nParams() ; ++m){
               systemPtr_->unitCell().SetParams( wCpArrays_[m] + lambda_ * dCpArrays_[m],m);
            }
            
            systemPtr_->unitCell().setLattice();
            systemPtr_->mixture().setupUnitCell(systemPtr_->unitCell());
            systemPtr_->basis().update(systemPtr_->unitCell());
            
            
            for(int m = 0; m < systemPtr_->unitCell().nParams() ; ++m){
               std::cout<<"Parameter"<<m<<"\t"<<"="
                        << std::setprecision (15)
                        <<systemPtr_->unitCell().params()[m]<<"\n";
            }

            gettimeofday(&timeEnd, &tz);       
            Log::file() <<"Parameter update time : "
                        <<Dbl( ( (double)(timeEnd.tv_sec - timeStart.tv_sec) + 
                                 (double)(timeEnd.tv_usec - timeStart.tv_usec) / 1000000.0 )
                               , 18, 11) <<'\n';

            

         }


      }
   }

}
}

#endif

