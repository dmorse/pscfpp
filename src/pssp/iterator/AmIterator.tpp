#ifndef PSSP_AM_ITERATOR_TPP
#define PSSP_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <pssp/System.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/containers/FArray.h>
#include <util/format/Dbl.h>
#include <ctime>
#include <cmath> 

namespace Pscf {
namespace Pssp 
{

   using namespace Util;

   #if 0
   template <int D>
   AmIterator<D>::AmIterator()
    : Iterator<D>(),
      epsilon_(0),
      lambda_(0),
      nHist_(0),
      maxHist_(0)
   { setClassName("AmIterator"); }
   #endif

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
      cell_ = 0; //default value (fixed cell)
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "domain", cell_); 
  }

   template <int D>
   void AmIterator<D>::allocate()
   {
      devHists_.allocate(maxHist_+1);
      omHists_.allocate(maxHist_+1);

      if (cell_){
         devCpHists_.allocate(maxHist_+1);
         CpHists_.allocate(maxHist_+1);
      }

      wArrays_.allocate(systemPtr_->mixture().nMonomer());
      dArrays_.allocate(systemPtr_->mixture().nMonomer());
      tempDev.allocate(systemPtr_->mixture().nMonomer());
     
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         wArrays_[i].allocate(systemPtr_->basis().nStar() - 1);
         dArrays_[i].allocate(systemPtr_->basis().nStar() - 1);
         tempDev[i].allocate(systemPtr_->basis().nStar() - 1);
      }


   }

   template <int D>
   int AmIterator<D>::solve()
   {
      //clock_t time_begin;
      //clock_t time_end;

      //solve the SCFT equations once
      //assumes basis.makeBasis() has been called
      //assumes AmIterator.allocate() has been called
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         systemPtr_->basis().convertFieldComponentsToDft(
                              systemPtr_->wField(i),
                              systemPtr_->wFieldDft(i));
         systemPtr_->fft().inverseTransform(systemPtr_->wFieldDft(i),
                                            systemPtr_->wFieldGrid(i));
      }

      systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                    systemPtr_->cFieldGrids());

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         systemPtr_->fft().forwardTransform(systemPtr_->cFieldGrid(i),
                                             systemPtr_->cFieldDft(i));
         systemPtr_->basis().convertFieldDftToComponents(
                              systemPtr_->cFieldDft(i),
                              systemPtr_->cField(i));
      }
       
      systemPtr_->computeFreeEnergy();
      Log::file()<<"Free Energy="<<systemPtr_->fHelmholtz();

      if (cell_){
         systemPtr_->mixture().computeStress();
         for (int m=0; m<(systemPtr_->unitCell()).nParameter() ; ++m){
               Log::file()<<"Stress"<<m<<"\t"<<"="<< systemPtr_->mixture().stress(m)<<"\n";
                 //Log::file()<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).params()[m]<<"\n";
               Log::file()<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).parameters()[m]<<"\n";
         }  
      } else {
         Log::file()<<std::endl;
      }

      //check for convergence else resolve SCFT equations with new Fields
      for (int itr = 1; itr <= maxItr_; ++itr) {
         
         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr-1;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         computeDeviation();

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<"  Iteration  "<<itr<<std::endl;

         if (isConverged()) {  
          Log::file()<<"----------CONVERGED----------"<< std::endl;
          if(!cell_){
             systemPtr_->mixture().computeStress();
          }
          for (int m=0; m<(systemPtr_->unitCell()).nParameter() ; ++m){
             Log::file()<<"Stress"<<m<<"\t"<<"="<< systemPtr_->mixture().stress(m)<<"\n";
             //Log::file()<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).params()[m]<<"\n";
             Log::file()<<"Parameter"<<m<<"\t"<<"="<<(systemPtr_->unitCell()).parameters()[m]<<"\n";
          }

              return 0;
     
         } else {
            if (itr <= maxHist_ + 1) {
               if (nHist_ > 0) {
                  invertMatrix_.allocate(nHist_, nHist_);
                  coeffs_.allocate(nHist_);
                  vM_.allocate(nHist_);
               }
            }

            minimizeCoeff(itr);
            buildOmega(itr);

            if (itr <= maxHist_) {
               //will deallocate when out of scope
               if (nHist_ > 0) {
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }

            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               systemPtr_->basis().convertFieldComponentsToDft( 
                                    systemPtr_->wField(j),
                                    systemPtr_->wFieldDft(j));
               systemPtr_->fft().inverseTransform(systemPtr_->wFieldDft(j), 
                                                   systemPtr_->wFieldGrid(j));
            }

            systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                          systemPtr_->cFieldGrids());

            if (cell_){
               systemPtr_->mixture().computeStress();
            }

            for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
               systemPtr_->fft().forwardTransform(systemPtr_->cFieldGrid(i),
                                                  systemPtr_->cFieldDft(i));
               systemPtr_->basis().convertFieldDftToComponents(
                                    systemPtr_->cFieldDft(i),
                                    systemPtr_->cField(i));
            }
         }

      }

      //should not reach here. iterated more than maxItr. Not converged
      return 1;
   }

   template <int D>
   void AmIterator<D>::computeDeviation()
   {

      omHists_.append(systemPtr_->wFields());

      if (cell_)
         //CpHists_.append((systemPtr_->unitCell()).params());
         CpHists_.append((systemPtr_->unitCell()).parameters());

      for (int i = 0 ; i < systemPtr_->mixture().nMonomer(); ++i) {
         for (int j = 0; j < systemPtr_->basis().nStar() - 1; ++j) {
            tempDev[i][j] = 0;
         }
      }

      DArray<double> temp;
      temp.allocate(systemPtr_->basis().nStar() - 1);

      #if 0
      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
         
         for (int j = 0; j < systemPtr_->basis().nStar() - 1; ++j) {
            temp[j] = 0;
         }

         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
               tempDev[i][k] += systemPtr_->interaction().chi(i,j) *
                              systemPtr_->cField(j)[k + 1];
               temp[k] += systemPtr_->wField(j)[k + 1];
            }
         }

         for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
            tempDev[i][k] += ((temp[k] / systemPtr_->mixture().nMonomer())
                             - systemPtr_->wField(i)[k + 1]);
         }
      }
      #endif

      for (int i = 0; i < systemPtr_->mixture().nMonomer(); ++i) {
    
         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
               tempDev[i][k] +=( (systemPtr_->interaction().chi(i,j) *
                              systemPtr_->cField(j)[k + 1]) - (systemPtr_->interaction().indemp(i,j) *
                              systemPtr_->wField(j)[k + 1]) ) ; 
            }   
         }   

      } 

      devHists_.append(tempDev);

      if (cell_){
         FArray<double, 6 > tempCp;
         for (int i = 0; i<(systemPtr_->unitCell()).nParameter() ; i++){
            tempCp [i] = -((systemPtr_->mixture()).stress(i));
         }
         devCpHists_.append(tempCp);
      }
   }
   
   template <int D>
   bool AmIterator<D>::isConverged()
   {
      double error;

     /// Error as defined in Matsen's Papers
      #if 0
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

      if (cell_){
         for ( int i = 0; i < (systemPtr_->unitCell()).nParameter() ; i++) {
            dError +=  devCpHists_[0][i] *  devCpHists_[0][i];
            //wError +=  (systemPtr_->unitCell()).params() [i] * (systemPtr_->unitCell()).params() [i];
            wError +=  (systemPtr_->unitCell()).parameters() [i] * (systemPtr_->unitCell()).parameters() [i];
         }
      }
      Log::file()<<" dError :"<<Dbl(dError)<<std::endl;
      Log::file()<<" wError :"<<Dbl(wError)<<std::endl;
      error = sqrt(dError / wError);
      #endif 

      /// Error by Max Residuals

      double temp1 = 0;
      double temp2 = 0;
      for ( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         for ( int j = 0; j < systemPtr_->basis().nStar() - 1; j++) {
            if (temp1 < fabs (devHists_[0][i][j]))
                temp1 = fabs (devHists_[0][i][j]);
         }   
      }
      Log::file()<<" SCF Error :"<<temp1<<std::endl;   
      error = temp1;

      if (cell_){
         for ( int i = 0; i < (systemPtr_->unitCell()).nParameter() ; i++) {
            if (temp2 < fabs (devCpHists_[0][i]))
                temp2 = fabs (devCpHists_[0][i]);
         }

         for (int m=0; m<(systemPtr_->unitCell()).nParameter() ; ++m){
            Log::file()<<" Stress "<<m<<" :"<< std::setprecision (15)<< systemPtr_->mixture().stress(m)<<"\n";
         }         
         error = (temp1>(100*temp2)) ? temp1 : (100*temp2);  
         // 100 is chose as stress rescale factor, seperate implementation of errors needs to be done
      }    

      // Depending on the value of cell_, residual errors?

      Log::file()<<" Error :"<<error<<std::endl;
      //Log::file()<<"---------------------"<<std::endl;
      if (error < epsilon_) {
         return true;
      } else {
         return false;
      }
   }

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
      //clock_t time_begin;
      //clock_t time_end;
      if (itr == 1) {
         //do nothing
      } else {

         double elm, elm_cp;

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

               if (cell_){
                  elm_cp = 0;
                  for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m){
                     elm_cp += ((devCpHists_[0][m] - devCpHists_[i+1][m]) * 
                                (devCpHists_[0][m] - devCpHists_[j+1][m])); 
                  }
                  invertMatrix_(i,j) += elm_cp;
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

            if (cell_){
               elm_cp = 0;
               for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m){
                  vM_[i] += ((devCpHists_[0][m] - devCpHists_[i+1][m]) *
                             (devCpHists_[0][m]));
               }
            }
         }

         if (itr == 2) {
            coeffs_[0] = vM_[0] / invertMatrix_(0,0);
         } else {
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);
            solver.solve(vM_, coeffs_);
         }       
      }
   }

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

         if (cell_){
            for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m){
 
                 parameters.append(CpHists_[0][m] +lambda_* devCpHists_[0][m]);
 
                 //(systemPtr_->unitCell()).SetParams( CpHists_[0][m] +lambda_* devCpHists_[0][m] ,m);
            }

                  (systemPtr_->unitCell()).setParameters(parameters);
                  (systemPtr_->unitCell()).setLattice();
                  systemPtr_->mixture().setupUnitCell(systemPtr_->unitCell());
                  systemPtr_->basis().update();

            for (int m=0; m<(systemPtr_->unitCell()).nParameter()  ; ++m){
               //Log::file()<<" Parameter "<<m<<" :"<<(systemPtr_->unitCell()).params()[m]<<"\n";
               Log::file()<<" Parameter "<<m<<" :"<<(systemPtr_->unitCell()).parameters()[m]<<"\n"; 
            }
         } 

      } else {
         for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
            for (int k = 0; k < systemPtr_->basis().nStar() - 1; ++k) {
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
            }
         }

         if(cell_){

            for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m){
               wCpArrays_[m] = CpHists_[0][m];
               dCpArrays_[m] = devCpHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m) {
                  wCpArrays_[m] += coeffs_[i] * ( CpHists_[i+1][m]-
                                                   CpHists_[0][m]);
                  dCpArrays_[m] += coeffs_[i] * ( devCpHists_[i+1][m]-
                                                   devCpHists_[0][m]);
               }
            } 
            for (int m = 0; m < (systemPtr_->unitCell()).nParameter() ; ++m){
               //(systemPtr_->unitCell()).SetParams( wCpArrays_[m] + lambda_ * dCpArrays_[m],m);

               parameters [m] = wCpArrays_[m] + lambda_ * dCpArrays_[m];
            }

            (systemPtr_->unitCell()).setParameters(parameters);
            (systemPtr_->unitCell()).setLattice();
            systemPtr_->mixture().setupUnitCell(systemPtr_->unitCell());
	    systemPtr_->basis().update();

            for (int m=0; m<(systemPtr_->unitCell()).nParameter() ; ++m){
	         // Log::file()<<" Parameter "<<m<<" :"<< std::setprecision (15)<<(systemPtr_->unitCell()).params()[m]<<"\n";
               Log::file()<<" Parameter "<<m<<" :"<< std::setprecision (15)<<(systemPtr_->unitCell()).parameters()[m]<<"\n";
            }


         }
      }
   }

}
}
#endif
