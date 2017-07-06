#ifndef PSSP_AM_ITERATOR_TPP
#define PSSP_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <util/format/Dbl.h>

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
   	/*for(int i = 0; i < nMonomer; i++){
   		monomerDevs_[i].allocate(nHist_, nGrid);
   	}*/

   	wArrays_.allocate(systemPtr_->mixture().nMonomer());
   	dArrays_.allocate(systemPtr_->mixture().nMonomer());

      for(int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         wArrays_[i].allocate(systemPtr_->basis().nStar());
         dArrays_[i].allocate(systemPtr_->basis().nStar());
      }

   }

   template <int D>
   int AmIterator<D>::solve()
   {
      //solve the SCFT equations once
      //assumes basis.makeBasis has been called
      //consider making dft arrays size of one monomer and reuse the same array
      for(int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         systemPtr_->basis().convertFieldComponentsToDft(
                              systemPtr_->wField(i),
                              systemPtr_->wFieldDft(i));
         systemPtr_->fft().inverseTransform( systemPtr_->wFieldDft(i),
                                             systemPtr_->wFieldGrid(i) );
      }

      systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                    systemPtr_->cFieldGrids());
      
      for(int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         systemPtr_->fft().forwardTransform( systemPtr_->cFieldGrid(i),
                                             systemPtr_->cFieldDft(i) );
         systemPtr_->basis().convertFieldDftToComponents(
                              systemPtr_->cFieldDft(i),
                              systemPtr_->cField(i) );
      }

      //check for convergence else resolve SCFT equations with new Fields
      for(int itr = 1; itr <= maxItr_; itr++) {
         
         if(itr <= maxHist_) {
            lambda_ = 1.0 - pow (0.9, itr);
            nHist_ = itr-1;
         }
         else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         computeDeviation();

         //std::cout<<"iterataion number: "<<itr<<std::endl;
         //std::cout<<"nHist_ : "<<nHist_<<std::endl;
         if(isConverged()) {
            return 0;
         }
         else {

            //resize history based matrix appropriately
            //consider making these working space local
            if(itr <= maxHist_ + 1 ) {
               if(nHist_ > 0){
                  invertMatrix_.allocate(nHist_, nHist_);
                  coeffs_.allocate(nHist_);
                  vM_.allocate(nHist_);
               }
            }
            minimizeCoeff(itr);
            buildOmega(itr);

            if(itr <= maxHist_){
               //will deallocate when out of scope
               if( nHist_ > 0){
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }

            for(int j = 0; j < systemPtr_->mixture().nMonomer(); j++) {
               systemPtr_->basis().convertFieldComponentsToDft( 
                                    systemPtr_->wField(j),
                                    systemPtr_->wFieldDft(j) );
               systemPtr_->fft().inverseTransform( systemPtr_->wFieldDft(j), 
                                                   systemPtr_->wFieldGrid(j) );
            }

            systemPtr_->mixture().compute(systemPtr_->wFieldGrids(), 
                                          systemPtr_->cFieldGrids());
         }

      }

      //should not reach here. iterated more than maxItr. Not converged
      return 1;
   }

   template <int D>
   void AmIterator<D>::computeDeviation()
	{

      omHists_.append(systemPtr_->wFields());

      DArray< DArray<double> > tempDev;
      tempDev.allocate(systemPtr_->mixture().nMonomer());

      for ( int i = 0 ; i < systemPtr_->mixture().nMonomer(); i++) {
         tempDev[i].allocate(systemPtr_->basis().nStar());

         for (int j = 0; j < systemPtr_->basis().nStar(); j++) {
            tempDev[i][j] = 0;
         }
      }
      //the form for this is slightly different for 3 species
      //Almost impossible to write good code here if using Interaction class
      //over ChiInteraction
      //I am not sure which one is faster but the first code uses cpu cache
      //to potentially speedup the code even with additional loops.
      DArray<double> temp;
      temp.allocate(systemPtr_->basis().nStar());

      for ( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         
         for(int j = 0; j < systemPtr_->basis().nStar(); j++) {
            temp[j] = 0;
         }

         for ( int j = 0; j < systemPtr_->mixture().nMonomer(); j++) {
            for ( int k = 0; k < systemPtr_->basis().nStar(); k++) {
               tempDev[i][k] += systemPtr_->interaction().chi(i,j) *
                              systemPtr_->cField(j)[k];
               temp[k] += systemPtr_->wField(j)[k];
            }
         }

         for ( int k = 0; k < systemPtr_->basis().nStar(); k++) {
            tempDev[i][k] += ( (temp[k] / systemPtr_->mixture().nMonomer())
                             - systemPtr_->wField(i)[k] );
         }
      }
      /*for ( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
         for (int j = 0; j < systemPtr_->basis().nStar(); j++) {
            temp = 0;
            for (int k = 0; k < systemPtr_->mixture().nMonomer(); k++) {
               tempDev[i][j] += (systemPtr_->interaction().chi(i,k) *
                                 systemPtr_->cField(k)[j] );
               temp += systemPtr_->wField(k)[j];
            }

            temp /= systemPtr_->mixture().nMonomer();

            tempDev[i][j] += temp - systemPtr_->wField(i)[j];
         }
      }*/
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
         for ( int j = 0; j < systemPtr_->basis().nStar(); j++) {
            //std::cout<<"nMonomer: "<<i<<" nStar : "<<j<<std::endl;
            dError += devHists_[0][i][j] * devHists_[0][i][j];
            //std::cout<<devHists_[0][i][j] * devHists_[0][i][j]<<std::endl;
            wError += systemPtr_->wField(i)[j] * systemPtr_->wField(i)[j];
            //std::cout<<systemPtr_->wField(i)[j] * systemPtr_->wField(i)[j]<<std::endl;
         }
      }
      error = sqrt(dError / wError);
      //std::cout<<"error "<<error<<std::endl;
      if( error < epsilon_)
         return true;
      else
         return false;
   }

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
      if(itr == 1){
         //do nothing
      }
      else {
         //std::cout<<"What is the size of hist + 1? "<<devHists_.size()<<std::endl;
         for(int i = 0; i < nHist_; i++) {
            vM_[i] = 0;
            for(int j = 0; j < nHist_; j++) {
               invertMatrix_(i,j) = 0;
            }
         }

         for(int i = 0; i < nHist_; i++) {

            for(int j = 0; j < nHist_; j++) {
               for(int k = 0; k < systemPtr_->mixture().nMonomer(); k++) {
                  for(int l = 0; l < systemPtr_->basis().nStar(); l++) {
                     invertMatrix_(i,j) += 
                        ( (devHists_[0][k][l] - devHists_[i+1][k][l]) *
                          (devHists_[0][k][l] - devHists_[j+1][k][l]) );
                  }
               }
            }

            for(int j = 0; j < systemPtr_->mixture().nMonomer(); j++) {
               for(int k = 0; k < systemPtr_->basis().nStar(); k++) {
                  vM_[i] += ( (devHists_[0][j][k] - devHists_[i+1][j][k]) *
                               devHists_[0][j][k] );
               }
            }
         }


         if(itr == 2){
            coeffs_[0] = vM_[0] / invertMatrix_(0,0);
         }
         else {
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(invertMatrix_);
            solver.solve(vM_, coeffs_);
            //std::cout<<"at least it's solved right?"<<std::endl;
         }       
      }
   }

   template <int D>
   void AmIterator<D>::buildOmega(int itr)
   {

      if(itr == 1){
         for( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
            for (int j = 0; j < systemPtr_->basis().nStar(); j++) {
               systemPtr_->wField(i)[j] = omHists_[0][i][j] + lambda_*devHists_[0][i][j];
            }
         }
      }
      else {
         //should be strictly correct. coeffs_ is a matrix of size 1 if itr ==2

         //std::cout<<"omHists_"<<std::endl;
         for( int j = 0; j < systemPtr_->mixture().nMonomer(); j++ ) {
            for( int k = 0; k < systemPtr_->basis().nStar(); k++) {
               wArrays_[j][k] = omHists_[0][j][k];
               dArrays_[j][k] = devHists_[0][j][k];
               //std::cout<<wArrays_[j][k]<<std::endl;
            }
         }

         //std::cout<<"wArrays and dArrays"<<std::endl;
         //std::cout<<"coeffs"<<std::endl;
         for(int i = 0; i < nHist_; i++) {
            for( int j = 0; j < systemPtr_->mixture().nMonomer(); j++) {
               for( int k = 0; k < systemPtr_->basis().nStar(); k++) {
                  wArrays_[j][k] += coeffs_[i] * ( omHists_[i+1][j][k] - 
                                                   omHists_[0][j][k] );
                  dArrays_[j][k] += coeffs_[i] * ( devHists_[i+1][j][k] - 
                                                   devHists_[0][j][k] );
                  //std::cout<<coeffs_[i]<<std::endl;
                  //std::cout<<wArrays_[j][k]<<std::endl;
                  //std::cout<<dArrays_[j][k]<<std::endl;
               }
            }
         }


         //std::cout<<"wFields"<<std::endl;
         for( int i = 0; i < systemPtr_->mixture().nMonomer(); i++) {
            for( int j = 0; j < systemPtr_->basis().nStar(); j++) {
              systemPtr_->wField(i)[j] = wArrays_[i][j] + lambda_ * dArrays_[i][j];
              //std::cout<<systemPtr_->wField(i)[j]<<std::endl;
            }
         }

      }
   }
}
}

#endif
