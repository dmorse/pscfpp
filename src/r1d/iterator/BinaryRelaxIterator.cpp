/*
 * PSCF - Polymer Self-Consistent Field Theory
 *
 * Copyright 2016 - 2022, The Regents of the University of Minnesota
 * Distributed under the terms of the GNU General Public License.
 */

#include "BinaryRelaxIterator.h"
#include <r1d/System.h>
#include <pscf/inter/Interaction.h>
#include <util/misc/Log.h>
#include <util/misc/Timer.h>
#include <math.h>

namespace Pscf {
namespace R1d
{
   using namespace Util;

   // Constructor
   BinaryRelaxIterator::BinaryRelaxIterator(System& system)
    : Iterator(system),
      epsilon_(0.0),
      lambdaPlus_(0.0),
      lambdaMinus_(0.0),
      maxItr_(200),
      isAllocated_(false),
      isCanonical_(true)
   {  setClassName("BinaryRelaxIterator"); }

   // Destructor
   BinaryRelaxIterator::~BinaryRelaxIterator()
   {}

   // Read parameter file block and allocate memory
   void BinaryRelaxIterator::readParameters(std::istream& in)
   {
      UTIL_CHECK(domain().nx() > 0);
      read(in, "epsilon", epsilon_);
      read(in, "maxIter", maxItr_);
      read(in, "lambdaPlus", lambdaPlus_);
      read(in, "lambdaMinus", lambdaMinus_);

      allocate();
   }

   // Allocate required memory (called by readParameters)
   void BinaryRelaxIterator::allocate()
   {
      if (isAllocated_) return;
      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);
      int nx = domain().nx();         // number of grid points
      UTIL_CHECK(nx > 0);
      cArray_.allocate(nm);
      wArray_.allocate(nm);
      dW_.allocate(nm);
      dWNew_.allocate(nm);
      wFieldsNew_.allocate(nm);
      cFieldsNew_.allocate(nm);
      for (int i = 0; i < nm; ++i) {
         wFieldsNew_[i].allocate(nx);
         cFieldsNew_[i].allocate(nx);
         dW_[i].allocate(nx);
         dWNew_[i].allocate(nx);
      }
      isAllocated_ = true;
   }
 
   void BinaryRelaxIterator::computeDW(Array<WField> const & wOld, 
                              Array<CField> const & cFields,
                              Array<WField> & dW,
                              double & dWNorm)
   
   {
      double dWm, dWp, c0, c1, w0, w1, wm;
      double dWmNorm = 0.0;
      double dWpNorm = 0.0;
      double chi = system().interaction().chi(0,1);
      int nx = domain().nx();        // number of grid points
      //For AB diblok
      for (int i = 0; i < nx; ++i) {
         c0 = cFields[0][i];
         c1 = cFields[1][i];
         w0 = wOld[0][i];
         w1 = wOld[1][i];
         wm = w1 - w0;
         dWp = c1 + c0 - 1.0;
         dWm = 0.5*( c0 - c1 - wm/chi);
         dWpNorm += dWp*dWp;
         dWmNorm += dWm*dWm;
         dWp *= lambdaPlus_;
         dWm *= lambdaMinus_;
         dW[0][i] = dWp -  dWm;
         dW[1][i] = dWp + dWm;
      }
      dWNorm = (dWpNorm + dWmNorm)/double(nx);
      dWNorm = sqrt(dWNorm);
   }

   void BinaryRelaxIterator::updateWFields(Array<WField> const & wOld,
                                  Array<WField> const & dW,
                                  Array<WField> & wNew)
   {
      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);
      int nx = domain().nx();        // number of grid points
      int i;                         // monomer index
      int j;                         // grid point index
      double w0, w1;
      
      // Add dW
      for (j = 0; j < nx; j++){
         w0 = wOld[0][j];
         w1 = wOld[1][j];
         wNew[0][j] = w0 + dW[0][j];
         wNew[1][j] = w1 + dW[1][j];
       }
    
      // If canonical, shift such that last element is exactly zero
      if (isCanonical_) {
         double shift = wNew[nm-1][nx-1];
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wNew[i][j] -= shift;
            }
         }
      }

   }

   int BinaryRelaxIterator::solve(bool isContinuation)
   {
      //Declare Timer
      Timer timerTotal;
      
      int nm = mixture().nMonomer();  // number of monomer types
      UTIL_CHECK(nm == 2);
      int np = mixture().nPolymer();  // number of polymer species
      int nx = domain().nx();         // number of grid points
      
      // Start overall timer 
      timerTotal.start();
      
      // Determine if isCanonical (iff all species ensembles are closed)
      isCanonical_ = true;
      Species::Ensemble ensemble;
      for (int i = 0; i < np; ++i) {
         ensemble = mixture().polymer(i).ensemble();
         if (ensemble == Species::Unknown) {
            UTIL_THROW("Unknown species ensemble");
         }
         if (ensemble == Species::Open) {
            isCanonical_ = false;
         }
      }
    
      // If isCanonical, shift so that last element is zero.
      // Note: This is one of the residuals in this case.
      if (isCanonical_) {
         double shift = wFields()[nm-1][nx-1];
         int i, j;
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wFields()[i][j] -= shift;
            }
         }
      }
      
      // Compute initial dWNorm.
      mixture().compute(system().wFields(), system().cFields());
      computeDW(system().wFields(), system().cFields(), dW_, dWNorm_);
    
      // Iterative loop
      int i, j, k;
      for (i = 0; i < maxItr_; ++i) {
         Log::file() << "iteration " << i
         << " , error = " << dWNorm_
         << std::endl;
        
         if (dWNorm_ < epsilon_) {
            // Stop timers
            timerTotal.stop();
            
            Log::file() << "The epsilon is " << epsilon_<< std::endl;
            Log::file() << "Converged" << std::endl;
            system().computeFreeEnergy();
            // Success
            Log::file() << "\n\n";
            // Output timing resultsl;
            Log::file() << "Total time:                             "  
                        << timerTotal.time()   << " s  "  << std::endl;
            Log::file() << "Average time cost of each iteration:    "  
                        << timerTotal.time()/i  << " s  " << std::endl;
            Log::file() << "\n\n";
            return 0;
         }

         // Try full Fd relaxation update
         updateWFields(system().wFields(), dW_, wFieldsNew_);
         mixture().compute(wFieldsNew_, cFieldsNew_);
         computeDW(wFieldsNew_, cFieldsNew_, dWNew_, dWNormNew_);
         
         // Decrease increment if necessary
         j = 0;
         while (dWNormNew_ > dWNorm_ && j < 3) {
            //double dWNormDecrease_;
            Log::file() << "      error = " 
                        << dWNormNew_ << ", decreasing increment" << std::endl;
            lambdaPlus_ *= 0.5;
            lambdaMinus_ *= 0.5;
            //Print lambdaPlus_ and lambdaMinus_ 
            Log::file() << "      lambdaPlus = " 
                        << lambdaPlus_ << std::endl;
            Log::file() << "      lambdaMinus = " 
                        << lambdaMinus_<< std::endl;                        
            computeDW(system().wFields(),system().cFields(), dWNew_, 
                      dWNormNew_);
            updateWFields(system().wFields(), dWNew_, wFieldsNew_);
            mixture().compute(wFieldsNew_, cFieldsNew_);
            ++j;
         }
         
         // Accept or reject update
         if (dWNormNew_ < dWNorm_) {
            // Update system fields
            for (j = 0; j < nm; ++j) {
               for (k = 0; k < nx; ++k) {
                  system().wField(j)[k] = wFieldsNew_[j][k];
                  system().cField(j)[k] = cFieldsNew_[j][k];
                  dW_[j][k] = dWNew_[j][k];
                }
            dWNorm_ = dWNormNew_;
            }
        } else {
            Log::file() << "Iteration failed, norm = "
            << dWNormNew_ << std::endl;
            break;
          }
        
      }
      // Failure
      return 1;
   }


} // namespace R1d
} // namespace Pscf
