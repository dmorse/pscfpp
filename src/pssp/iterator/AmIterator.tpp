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
#include <util/format/Dbl.h>
#include <ctime>

namespace Pscf {
namespace Pssp 
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
      //consider making dft arrays size of one monomer and reuse the same array
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

      //check for convergence else resolve SCFT equations with new Fields
      for (int itr = 1; itr <= maxItr_; ++itr) {
         
         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr-1;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         //time_begin = clock();
         computeDeviation();
         //time_end = clock();
         //std::cout<<" Time for computeDeviation ="
         //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;

         std::cout<<"  Iteration  "<<itr<<std::endl;
         if (isConverged()) {
            return 0;
         } else {
            //resize history based matrix appropriately
            //consider making these working space local
            if (itr <= maxHist_ + 1) {
               if (nHist_ > 0) {
                  invertMatrix_.allocate(nHist_, nHist_);
                  coeffs_.allocate(nHist_);
                  vM_.allocate(nHist_);
               }
            }
            //time_begin = clock();
            minimizeCoeff(itr);
            //time_end = clock();
            //std::cout<<" Time for minimizeCoeff ="
            //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;
            //time_begin = clock();
            buildOmega(itr);
            //time_end = clock();
            //std::cout<<" Time for buildOmega ="
            //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;

            if (itr <= maxHist_) {
               //will deallocate when out of scope
               if (nHist_ > 0) {
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }

            //time_begin = clock();
            for (int j = 0; j < systemPtr_->mixture().nMonomer(); ++j) {
               systemPtr_->basis().convertFieldComponentsToDft( 
                                    systemPtr_->wField(j),
                                    systemPtr_->wFieldDft(j));
               systemPtr_->fft().inverseTransform(systemPtr_->wFieldDft(j), 
                                                   systemPtr_->wFieldGrid(j));
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
            //time_end = clock();
            //std::cout<<" Time for within loop transform ="
            //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;
         }

      }

      //should not reach here. iterated more than maxItr. Not converged
      return 1;
   }

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
      std::cout<<" dError :"<<Dbl(dError)<<std::endl;
      std::cout<<" wError :"<<Dbl(wError)<<std::endl;
      error = sqrt(dError / wError);
      std::cout<<"  Error  :"<<Dbl(error)<<std::endl;
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

         /*for (int i = 0; i < nHist_; ++i) {
            for (int j = i; j < nHist_; ++j) {


               invertMatrix_(i,j) = 0;

               for (int k = 0; k < systemPtr_->mixture().nMonomer(); ++k) {
                  for (int l = 0; l < systemPtr_->basis().nStar() - 1; ++l) {
                     invertMatrix_(i,j) +=
                        (  (devHists_[0][k][l] - devHists_[i+1][k][l]) *
                           (devHists_[0][k][l] - devHists_[j+1][k][l]) );
                  }
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
         }*/

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
            //time_begin = clock();
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);
            solver.solve(vM_, coeffs_);
            //time_end = clock();
            //std::cout<<" nHist_ is "<<nHist_<<std::endl;
            //std::cout<<" Time for LUSolver ="
            //<< Dbl((float)(time_end - time_begin)/CLOCKS_PER_SEC,18,11)<<std::endl;
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
            }
         }

      }
   }
}
}

#endif
