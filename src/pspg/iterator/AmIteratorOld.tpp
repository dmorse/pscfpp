#ifndef PSPG_AM_ITERATOR_OLD_TPP
#define PSPG_AM_ITERATOR_OLD_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorOld.h"
#include <pspg/System.h>
#include <pspg/math/GpuResources.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/format/Dbl.h>
#include <util/containers/FArray.h>
#include <util/misc/Timer.h>
#include <sys/time.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   template <int D>
   AmIteratorOld<D>::AmIteratorOld(System<D>& system)
      : Iterator<D>(system),
        epsilon_(0),
        lambda_(0),
        nHist_(0),
        maxHist_(0),
        isFlexible_(0)
   { setClassName("AmIteratorOld"); }

   template <int D>
   AmIteratorOld<D>::~AmIteratorOld()
   {
      delete[] temp_;
      cudaFree(d_temp_);
   }

   template <int D>
   void AmIteratorOld<D>::readParameters(std::istream& in)
   {   
      isFlexible_ = 0;
      errorType_ = "normResid"; // default type of error
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "errorType", errorType_);

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);

      if (!(errorType_ == "normResid" || errorType_ == "maxResid")) {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }
   }

   template <int D>
   void AmIteratorOld<D>::setup()
   {
      // GPU resources
      const int size = system().mesh().size();
      int NUMBER_OF_BLOCKS, THREADS_PER_BLOCK;
      ThreadGrid::setThreadsLogical(size, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

      d_resHists_.allocate(maxHist_ + 1);
      d_omHists_.allocate(maxHist_ + 1);

      if (isFlexible_) {
         devCpHists_.allocate(maxHist_+1);
         CpHists_.allocate(maxHist_+1);
      }

      wArrays_.allocate(system().mixture().nMonomer());
      dArrays_.allocate(system().mixture().nMonomer());
      tempDev.allocate(system().mixture().nMonomer());

      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
          wArrays_[i].allocate(system().mesh().size());
          dArrays_[i].allocate(system().mesh().size());
          tempDev[i].allocate(system().mesh().size());
      }
      
      histMat_.allocate(maxHist_ + 1);
      //allocate d_temp_ here i suppose
      cudaMalloc((void**)&d_temp_, NUMBER_OF_BLOCKS * sizeof(cudaReal));
      temp_ = new cudaReal[NUMBER_OF_BLOCKS];
   }

   template <int D>
   int AmIteratorOld<D>::solve()
   {
      
      // Define Timer objects
      Timer solverTimer;
      Timer stressTimer;
      Timer updateTimer;
      Timer::TimePoint now;
      bool done;

      // Solve MDE for initial state
      solverTimer.start();
      system().mixture().compute(system().wFieldsRGrid(),
         system().cFieldsRGrid());
      now = Timer::now();
      solverTimer.stop(now);

      // Compute stress for initial state
      if (isFlexible_) {
         stressTimer.start(now);
         system().mixture().computeStress(system().wavelist());
         for (int m = 0; m < system().unitCell().nParameter() ; ++m){
            Log::file() << "Stress    " << m << " = "
                        << system().mixture().stress(m)<<"\n";
         }
         for (int m = 0; m < system().unitCell().nParameter() ; ++m){
            Log::file() << "Parameter " << m << " = "
                        << (system().unitCell()).parameter(m)<<"\n";
         }
         now = Timer::now();
         stressTimer.stop(now);
      }

      // Anderson-Mixing iterative loop
      int itr;
      for (itr = 1; itr <= maxItr_; ++itr) {
         updateTimer.start(now);
        
         Log::file() << "---------------------" << std::endl;
         Log::file() << " Iteration  " << itr << std::endl;
         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr - 1;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         computeDeviation();
         done = isConverged();

         if (done) {
            updateTimer.stop();

            if (itr > maxHist_ + 1) {
               invertMatrix_.deallocate();
               coeffs_.deallocate();
               vM_.deallocate();
            }

            Log::file() << "------- CONVERGED ---------"<< std::endl;

            // Output final timing results
            double updateTime = updateTimer.time();
            double solverTime = solverTimer.time();
            double stressTime = 0.0;
            double totalTime = updateTime + solverTime;
            if (isFlexible_) {
               stressTime = stressTimer.time();
               totalTime += stressTime;
            }
            Log::file() << "\n";
            Log::file() << "Iterator times contributions:\n";
            Log::file() << "\n";
            Log::file() << "solver time  = " << solverTime  << " s,  "
                        << solverTime/totalTime << "\n";
            Log::file() << "stress time  = " << stressTime  << " s,  "
                        << stressTime/totalTime << "\n";
            Log::file() << "update time  = "  << updateTime  << " s,  "
                        << updateTime/totalTime << "\n";
            Log::file() << "total time   = "  << totalTime   << " s  ";
            Log::file() << "\n\n";
           
            if (isFlexible_) {
               Log::file() << "\n";
               Log::file() << "Final stress values:" << "\n";
               for (int m = 0; m < system().unitCell().nParameter() ; ++m){
                  Log::file() << "Stress    " << m << " = "
                              << system().mixture().stress(m)<<"\n";
               }
               Log::file() << "\n";
               Log::file() << "Final unit cell parameter values:" << "\n";
               for (int m = 0; m < system().unitCell().nParameter() ; ++m){
                  Log::file() << "Parameter " << m << " = "
                              << (system().unitCell()).parameter(m)<<"\n";
               }
               Log::file() << "\n";
            }
            return 0;

         } else {

            // Resize history based matrix appropriately
            // consider making these working space local

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
            now = Timer::now();
            updateTimer.stop(now);

            // Solve MDE
            solverTimer.start(now);
            system().mixture().compute(system().wFieldsRGrid(),
                                          system().cFieldsRGrid());
            now = Timer::now();
            solverTimer.stop(now);
     
            if (isFlexible_) {
               stressTimer.start(now);
               system().mixture().computeStress(system().wavelist());
               for (int m = 0; m < system().unitCell().nParameter() ; ++m){
                  Log::file() << "Stress    " << m << " = "
                              << system().mixture().stress(m)<<"\n";
               }
               for (int m = 0; m < system().unitCell().nParameter() ; ++m){
                  Log::file() << "Parameter " << m << " = "
                              << (system().unitCell()).parameter(m)<<"\n";
               }
               now = Timer::now();
               stressTimer.stop(now);
            }
            
         }
         
      }

      if (itr > maxHist_ + 1) {
         invertMatrix_.deallocate();
         coeffs_.deallocate();
         vM_.deallocate();
      }
 
      // Failure: Not converged after maxItr iterations.
      return 1;
   }

   template <int D>
   void AmIteratorOld<D>::computeDeviation()
   {
      // GPU resources
      const int size = system().mesh().size();
      int NUMBER_OF_BLOCKS, THREADS_PER_BLOCK;
      ThreadGrid::setThreadsLogical(size, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

      //need to average
      float average = 0;
      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         average += gpuSum(system().wFieldRGrid(i).cDField(), size);
      }
      average /= (system().mixture().nMonomer() * size);
      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         subtractUniform <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (system().wFieldRGrid(i).cDField(), average, size);
      }

      d_omHists_.append(system().wFieldsRGrid());

      if (isFlexible_) {
         CpHists_.append(system().unitCell().parameters());
      }

      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         assignUniformReal <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (tempDev[i].cDField(), 0, size);
      }
      
      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         for (int j = 0; j < system().mixture().nMonomer(); ++j) {
            pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> 
               (tempDev[i].cDField(),
                system().cFieldRGrid(j).cDField(),
                system().interaction().chi(i, j),
                size);
            //this is a good add but i dont necessarily know if its right
            pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (tempDev[i].cDField(),
                                                                            system().wFieldRGrid(j).cDField(),
                                                                            -system().interaction().idemp(i, j),
                                                                            size); 

         }
         
      }

      float sum_chi_inv = (float) system().interaction().sum_inv();

      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         
         pointWiseSubtractFloat <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (tempDev[i].cDField(),
                                                                             1/sum_chi_inv, 
                                                                             size);
      }

      d_resHists_.append(tempDev);

      if (isFlexible_) {
         FArray<double, 6> tempCp;
         for (int i = 0; i<(system().unitCell()).nParameter(); i++){
            //format????
            tempCp [i] = -((system().mixture()).stress(i));
         }
         devCpHists_.append(tempCp);
      }
   }

   template <int D>
   bool AmIteratorOld<D>::isConverged()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int size = system().mesh().size();

      double matsenError;
      double dError = 0;
      double wError = 0;
      
      for (int i = 0; i < nMonomer; ++i) {
         dError += gpuInnerProduct(d_resHists_[0][i].cDField(), d_resHists_[0][i].cDField(), size);
         wError += gpuInnerProduct(system().wFieldRGrid(i).cDField(), system().wFieldRGrid(i).cDField(), size);
      }

      if (isFlexible_) {
         for ( int i = 0; i < nParameter; i++) {
            dError +=  devCpHists_[0][i] *  devCpHists_[0][i];
            wError +=  system().unitCell().parameter(i) * system().unitCell().parameter(i);
         }
      }

      Log::file() << " dError :" << Dbl(dError) << '\n';
      Log::file() << " wError :" << Dbl(wError) << '\n';
      matsenError = sqrt(dError / wError);
      Log::file() << "  Error :" << Dbl(matsenError) << '\n';

      // Find max residuals. 
      cudaReal* tempRes = new cudaReal[nMonomer*size + nParameter];
      double maxSCF=0.0, maxStress=0.0, maxRes=0.0;
      for (int i = 0; i < nMonomer; i++ ) {
         cudaMemcpy(&tempRes[i*size], d_resHists_[0][i].cDField(), size*sizeof(cudaReal), cudaMemcpyDeviceToHost);
      }

      for (int i = 0; i < nMonomer*size; i++ ) {
         if (fabs(tempRes[i] > maxSCF)) maxSCF = tempRes[i]; 
      }
      maxRes = maxSCF;
      Log::file() << "  Max SCF Residual :" << Dbl(maxRes) << '\n';

      if (isFlexible_) {
         for (int i = 0; i < nParameter; i++ ) {
            if (fabs(devCpHists_[0][i]) > maxStress) maxStress = devCpHists_[0][i];
         }
         // compare SCF and stress residuals 
         if (maxStress*10 > maxRes) maxRes = maxStress;
         Log::file() << "  Max Stress Residual :" << Dbl(maxStress*10) << '\n';
      }
            
      // return convergence boolean for chosen error type
      if (errorType_ == "normResid") {
         return matsenError < epsilon_;
      } else if (errorType_ == "maxResid") {
         return maxRes < epsilon_;
      } else {
         UTIL_THROW("Invalid iterator error type in parameter file.");
         return false;
      }

            
   }

   template <int D>
   int AmIteratorOld<D>::minimizeCoeff(int itr)
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
               for (int j = 0; j < system().mixture().nMonomer(); ++j) {
                  elm += gpuInnerProduct(d_resHists_[0][j].cDField(), d_resHists_[i][j].cDField(), system().mesh().size());
               }
               histMat_.evaluate(elm, nHist_, i);
            }
         }

         //build new Umn array
         for (int i = 0; i < nHist_; ++i) {
            for (int j = i; j < nHist_; ++j) {
               invertMatrix_(i, j) = histMat_.makeUmn(i, j, nHist_);

               if (isFlexible_) {
                  elm_cp = 0;
                  for (int m = 0; m < system().unitCell().nParameter(); ++m){
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

            if (isFlexible_) {
               elm_cp = 0;
               for (int m = 0; m < system().unitCell().nParameter(); ++m){
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
   void AmIteratorOld<D>::buildOmega(int itr)
   {

      // GPU resources
      const int size = system().mesh().size();
      int NUMBER_OF_BLOCKS, THREADS_PER_BLOCK;
      ThreadGrid::setThreadsLogical(size, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

      if (itr == 1) {
         for (int i = 0; i < system().mixture().nMonomer(); ++i) {
            assignReal <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (system().wFieldRGrid(i).cDField(),
            d_omHists_[0][i].cDField(), size);
            pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (system().wFieldRGrid(i).cDField(),
            d_resHists_[0][i].cDField(), lambda_, size);
         }

         if (isFlexible_) {
            cellParameters_.clear();
            for (int m = 0; m < (system().unitCell()).nParameter() ; ++m){
               cellParameters_.append(CpHists_[0][m] +lambda_* devCpHists_[0][m]);
            }
            system().unitCell().setParameters(cellParameters_);            
            system().mixture().setupUnitCell(system().unitCell(), system().wavelist());
            system().wavelist().computedKSq(system().unitCell());
         }

      } else {
         //should be strictly correct. coeffs_ is a vector of size 1 if itr ==2

         for (int j = 0; j < system().mixture().nMonomer(); ++j) {
            assignReal <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (wArrays_[j].cDField(),
               d_omHists_[0][j].cDField(), size);
            assignReal <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (dArrays_[j].cDField(),
               d_resHists_[0][j].cDField(), size);
         }

         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < system().mixture().nMonomer(); ++j) {
               //wArrays
               pointWiseBinarySubtract <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (d_omHists_[i + 1][j].cDField(),
                  d_omHists_[0][j].cDField(), tempDev[0].cDField(),
                  size);
               pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (wArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i], size);

               //dArrays
               pointWiseBinarySubtract <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (d_resHists_[i + 1][j].cDField(),
                  d_resHists_[0][j].cDField(), tempDev[0].cDField(),
                  size);
               pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (dArrays_[j].cDField(),
                  tempDev[0].cDField(), coeffs_[i], size);
            }
         }
         
         for (int i = 0; i < system().mixture().nMonomer(); ++i) {
            assignReal <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (system().wFieldRGrid(i).cDField(),
               wArrays_[i].cDField(), size);
            pointWiseAddScale <<< NUMBER_OF_BLOCKS, THREADS_PER_BLOCK >>> (system().wFieldRGrid(i).cDField(),
               dArrays_[i].cDField(), lambda_, size);
         }

         if (isFlexible_) {
            
            for (int m = 0; m < system().unitCell().nParameter() ; ++m){
               wCpArrays_[m] = CpHists_[0][m];
               dCpArrays_[m] = devCpHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < system().unitCell().nParameter() ; ++m) {
                  wCpArrays_[m] += coeffs_[i] * ( CpHists_[i+1][m]-
                                                  CpHists_[0][m]);
                  dCpArrays_[m] += coeffs_[i] * ( devCpHists_[i+1][m]-
                                                  devCpHists_[0][m]);
               }
            } 

            cellParameters_.clear();
            for (int m = 0; m < system().unitCell().nParameter() ; ++m){               
               cellParameters_.append(wCpArrays_[m] + lambda_* dCpArrays_[m]);
            }
            
            system().unitCell().setParameters(cellParameters_);            
            system().mixture().setupUnitCell(system().unitCell(), system().wavelist());
            system().wavelist().computedKSq(system().unitCell());
            
         }


      }
   }

}
}
#endif
