#ifndef PSPC_AM_ITERATOR_TPP
#define PSPC_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <pspc/System.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/containers/FArray.h>
#include <util/format/Dbl.h>
#include <util/misc/Timer.h>
#include <cmath>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AmIterator<D>::AmIterator(System<D>& system)
    : Iterator<D>(system),
      epsilon_(0),
      lambda_(0),
      nHist_(0),
      maxHist_(0)
   {  setClassName("AmIterator"); }

   /*
   * Destructor
   */
   template <int D>
   AmIterator<D>::~AmIterator()
   {}

   /*
   * Read parameter file block.
   */
   template <int D>
   void AmIterator<D>::readParameters(std::istream& in)
   {
      isFlexible_ = 0; // default value (fixed cell)
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "isFlexible", isFlexible_);
   }

   /*
   * Setup and allocate memory required by iterator.
   */
   template <int D>
   void AmIterator<D>::setup()
   {
      // Allocate ring buffers
      resHists_.allocate(maxHist_+1);
      wHists_.allocate(maxHist_+1);

      if (isFlexible_) {
         stressHists_.allocate(maxHist_+1);
         cellParamHists_.allocate(maxHist_+1);
      }

      // Allocate outer arrays used in iteration (number of monomers)
      int nMonomer = system().mixture().nMonomer();
      wArrays_.allocate(nMonomer);
      dArrays_.allocate(nMonomer);
      resArrays_.allocate(nMonomer);

      // Determine number of basis functions (nBasis - 1 if canonical)
      if (isCanonical()) {
         shift_ = 1;
      } else {
         shift_ = 0;
      }
      
      // Allocate inner arrays with number of basis functions
      int nBasis = system().basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         wArrays_[i].allocate(nBasis);
         dArrays_[i].allocate(nBasis);
         resArrays_[i].allocate(nBasis);
      }
   }

   /*
   * Solve iteratively.
   */
   template <int D>
   int AmIterator<D>::solve()
   {
      // Preconditions:
      UTIL_CHECK(system().hasWFields());
      // Assumes basis.makeBasis() has been called
      // Assumes AmIterator.allocate() has been called
      // TODO: Check these conditions on entry

      Timer convertTimer;
      Timer solverTimer;
      Timer stressTimer;
      Timer updateTimer;
      Timer::TimePoint now;
      bool done;

      // Solve MDE for initial state
      solverTimer.start();
      system().compute();
      now = Timer::now();
      solverTimer.stop(now);

      // Compute initial stress if needed
      if (isFlexible_) {
         stressTimer.start(now);
         system().mixture().computeStress();
         now = Timer::now();
         stressTimer.stop(now);
      }

      // Iterative loop
      for (int itr = 1; itr <= maxItr_; ++itr) {

         updateTimer.start(now);

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<" Iteration  "<<itr<<std::endl;

         if (itr <= maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr);
            nHist_ = itr-1;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }
         computeResidual();

         // Test for convergence
         done = isConverged();

         if (done) {

            updateTimer.stop();
            Log::file() << "----------CONVERGED----------"<< std::endl;

            // Output timing results
            double solverTime = solverTimer.time();
            double convertTime = convertTimer.time();
            double updateTime = updateTimer.time();
            double totalTime = updateTime + convertTime + solverTime;
            double stressTime = 0.0;
            if (isFlexible_) {
               stressTime = stressTimer.time();
               totalTime += stressTime;
            }
            Log::file() << "\n";
            Log::file() << "Iterator times contributions:\n";
            Log::file() << "\n";
            Log::file() << "solver time  = " << solverTime  << " s,  "
                        << solverTime/totalTime << "\n";
            if (isFlexible_) {
               Log::file() << "stress time  = " << stressTime  << " s,  "
                           << stressTime/totalTime << "\n";
            }
            Log::file() << "convert time = " << convertTime << " s,  "
                        << convertTime/totalTime << "\n";
            Log::file() << "update time  = "  << updateTime  << " s,  "
                        << updateTime/totalTime << "\n";
            Log::file() << "total time   = "  << totalTime   << " s  ";
            Log::file() << "\n\n";

            // If the unit cell is rigid, compute and output final stress 
            if (!isFlexible_) {
               system().mixture().computeStress();
               Log::file() << "Final stress:" << "\n";
               for (int m=0; m < system().unitCell().nParameter(); ++m) {
                  Log::file() << "Stress  "<< m << "   = "
                              << Dbl(system().mixture().stress(m)) 
                              << "\n";
               }
               Log::file() << "\n";
            }

            // Successful completion (i.e., converged within tolerance)
            cleanUp();
            return 0;

         } else {
            if (itr <= maxHist_ + 1) {
               if (nHist_ > 0) {
                  U_.allocate(nHist_, nHist_);
                  coeffs_.allocate(nHist_);
                  v_.allocate(nHist_);
               }
            }
            minimizeCoeff(itr);
            buildOmega(itr);

            if (itr <= maxHist_) {
               if (nHist_ > 0) {
                  U_.deallocate();
                  coeffs_.deallocate();
                  v_.deallocate();
               }
            }
            now = Timer::now();
            updateTimer.stop(now);

            // Solve MDE
            solverTimer.start(now);
            system().compute();
            now = Timer::now();
            solverTimer.stop(now);

            // Compute stress if needed
            if (isFlexible_) {
               stressTimer.start(now);
               system().mixture().computeStress();
               now = Timer::now();
               stressTimer.stop(now);
            }

         }

      }
      // Failure: iteration counter itr reached maxItr without converging
      cleanUp();
      return 1;
   }

   template <int D>
   void AmIterator<D>::computeResidual()
   {
      // Relevant quantities
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int nBasis = system().basis().nBasis();

      // Store current w field in history ringbuffer 
      wHists_.append(system().wFields());
      
      // If variable unit cell, store current unit cell parameters 
      if (isFlexible_) {
         cellParamHists_.append(system().unitCell().parameters());
      }

      // Initialize temporary residuals workspace 
      for (int i = 0 ; i < nMonomer; ++i) {
         for (int k = 0; k < nBasis; ++k) {
            resArrays_[i][k] = 0.0;
         }
      }
      
      // Compute residual vector
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = shift_; k < nBasis; ++k) {
               resArrays_[i][k] +=
                 ( system().interaction().chi(i,j)*system().cField(j)[k]
                 - system().interaction().idemp(i,j)*system().wField(j)[k]);
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (shift_ == 0) {
         for (int i = 0; i < nMonomer; ++i) {
            resArrays_[i][0] -= 1.0/system().interaction().sum_inv();
         }
      }

      // Store residuals
      resHists_.append(resArrays_);

      // If variable unit cell, store stress
      if (isFlexible_) {
         FArray<double, 6 > tempCp;
         for (int i = 0; i < nParameter ; i++) {
            tempCp [i] = -((system().mixture()).stress(i));
         }
         stressHists_.append(tempCp);
      }
   }

   template <int D>
   bool AmIterator<D>::isConverged()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int nBasis = system().basis().nBasis();

      double error = 0.0;

      #if 0
      // Error as defined in Matsen's Papers
      double dError = 0;
      double wError = 0;
      for ( int i = 0; i < nMonomer; i++) {
         for ( int j = shift_; j < nBasis; j++) {
            dError += resHists_[0][i][j] * resHists_[0][i][j];
            wError += system().wField(i)[j] * system().wField(i)[j];
         }
      }

      if (isFlexible_) {
         for ( int i = 0; i < nParameter ; i++) {
            dError += stressHists_[0][i] *  stressHists_[0][i];
            wError += system().unitCell().parameters()[i]
                     *system().unitCell().parameters()[i];
         }
      }
      Log::file() << " dError :" << Dbl(dError)<<std::endl;
      Log::file() << " wError :" << Dbl(wError)<<std::endl;
      error = sqrt(dError / wError);
      #endif

      // Error by max SCF residual
      double errSCF = 0.0;
      for ( int i = 0; i < nMonomer; i++) {
         for ( int j = shift_; j < nBasis; j++) {
            if (errSCF < fabs (resHists_[0][i][j]))
                errSCF = fabs (resHists_[0][i][j]);
         }
      }
      Log::file() << "SCF Error   = " << Dbl(errSCF) << std::endl;
      error = errSCF;

      // Error by max Stress residual
      if (isFlexible_) {
         double errStress = 0.0;
         for ( int i = 0; i < nParameter ; i++) {
            if (errStress < fabs (stressHists_[0][i])) {
                errStress = fabs (stressHists_[0][i]);
            }
         }
         // Output current stress values
         for (int m=0;  m < nParameter ; ++m) {
            Log::file() << "Stress  "<< m << "   = "
                        << Dbl(system().mixture().stress(m)) <<"\n";
         }
         // Determine if stress error, scaled by a factor, is larger
         // than the SCF error. If so, use it.  
         double scaleStress = 100.0;
         if (error < scaleStress*errStress) {
            error = scaleStress*errStress;
         }
      }
      Log::file() << "Error       = " << Dbl(error) << std::endl;

      // Output current unit cell parameter values
      if (isFlexible_) {
         for (int m=0; m < nParameter ; ++m) {
               Log::file() << "Parameter " << m << " = "
                           << Dbl(system().unitCell().parameters()[m])
                           << "\n";
         }
      }

      #if 0
      Log::file() << "n=0 Components of w basis" << std::endl;
      for (int i=0; i < nMonomer; ++i) {
         Log::file() << system().wField(i)[0] << std::endl;
      }
      #endif

      // Check if total error is below tolerance
      return error < epsilon_;
   }

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int nBasis = system().basis().nBasis();

      if (itr == 1) {
         //do nothing
      } else {         
         // Compute U matrix, as described in Arora 2017.
         for (int m = 0; m < nHist_; ++m) {
            for (int n = m; n < nHist_; ++n) {
               // Initialize U element value
               U_(m,n) = 0;
               // Compute dot products of differences of residual vectors 
               for (int i = 0; i < nMonomer; ++i) {
                  for (int k = shift_; k < nBasis; ++k) {
                     U_(m,n) +=
                            ((resHists_[0][i][k] - resHists_[m+1][i][k])*
                             (resHists_[0][i][k] - resHists_[n+1][i][k]));
                  }
               }
               // Include the stress residual contribution, if flexible cell
               if (isFlexible_) {
                  for (int p = 0; p < nParameter ; ++p) {
                     U_(m,n) += ((stressHists_[0][p] - stressHists_[m+1][p])*
                                (stressHists_[0][p] - stressHists_[n+1][p]));
                  }
               }
               U_(n,m) = U_(m,n);
            }
         }

         // Compute v vector, as described in Arora 2017.
         for (int m = 0; m < nHist_; ++m) {
            // Initialize v element value. 
            v_[m] = 0;
            // dot product of residual vectors
            for (int i = 0; i < nMonomer; ++i) {
               for (int k = shift_; k < nBasis; ++k) {
                  v_[m] += ( (resHists_[0][i][k] - resHists_[m+1][i][k]) *
                               resHists_[0][i][k] );
               }
            }
            // contributions of stress residuals if flexible cell
            if (isFlexible_) {
               for (int p = 0; p < nParameter ; ++p) {
                  v_[m] += ((stressHists_[0][p] - stressHists_[m+1][p]) *
                             (stressHists_[0][p]));
               }
            }
         }
         
         if (itr == 2) {
            // solve explicitly for coefficient
            coeffs_[0] = v_[0] / U_(0,0);
         } else {
            // numerically solve for coefficients
            LuSolver solver;
            solver.allocate(nHist_);
            solver.computeLU(U_);
            solver.solve(v_, coeffs_);
         }

         // output U for fun
         std::cout << "\n";
         for (int m = 0; m < nHist_; ++m) {
            for (int n = 0; n < nHist_; ++n) {
               std::cout << U_(m,n) << "  ";
            }
            std::cout << "\n";
         }

         std::cout << "\n";
         for (int m = 0; m < nHist_; ++m) {
            std::cout << v_[m] << std::endl;
         }
         std::cout << "\n";

         for (int m = 0; m < nHist_; ++m) {
            std::cout << coeffs_[m] << std::endl;
         }
         std::cout << "\n";
      }
   }

   template <int D>
   void AmIterator<D>::buildOmega(int itr)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int nBasis = system().basis().nBasis();

      if (itr == 1) { // if only 1 historical solutions

         // Update omega field with SCF residuals
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = shift_; k < nBasis; ++k) {
               wArrays_[i][k]
                      = wHists_[0][i][k] + lambda_*resHists_[0][i][k];
            }
         }
         // If canonical, set the homogeneous components explicitly
         if (shift_ == 1) {
            for (int i = 0; i < nMonomer; ++i) {
               dArrays_[i][0] = 0.0;
               wArrays_[i][0] = 0.0;
               for (int j = 0; j < nMonomer; ++j) {
                  wArrays_[i][0] += 
                     system().interaction().chi(i,j)*system().cField(j)[0];
               }
            }
         }
         system().setWBasis(wArrays_);

         // Update unit cell parameters with stress residuals
         if (isFlexible_) {
            parameters_.clear();
            for (int m = 0; m < nParameter ; ++m) {
               parameters_.append(cellParamHists_[0][m]
                              + lambda_* stressHists_[0][m]);

            }
            system().setUnitCell(parameters_);
         }

      } else { // if at least 2 historical results
         // contribution of the last solution
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = shift_; k < nBasis; ++k) {
               wArrays_[j][k] = wHists_[0][j][k];
               dArrays_[j][k] = resHists_[0][j][k];
            }
         }
         // mixing in historical solutions
         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < nMonomer; ++j) {
               for (int k = shift_; k < nBasis; ++k) {
                  wArrays_[j][k] += coeffs_[i] * ( wHists_[i+1][j][k] -
                                                   wHists_[0][j][k] );
                  dArrays_[j][k] += coeffs_[i] * ( resHists_[i+1][j][k] -
                                                   resHists_[0][j][k] );
               }
            }
         }
         // If isCanonical, set the homogeneous components explicitly
         if (shift_ == 1) {
            for (int i = 0; i < nMonomer; ++i) {
               dArrays_[i][0] = 0.0;
               wArrays_[i][0] = 0.0;
               for (int j = 0; j < nMonomer; ++j) {
                  wArrays_[i][0] += 
                     system().interaction().chi(i,j)*system().cField(j)[0];
               }
            }
         }
         // Adding in the estimated error of the predicted field
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               wArrays_[i][k] += lambda_ * dArrays_[i][k];
            }
         }
         // Updating the system
         system().setWBasis(wArrays_);

         std::cout << "Updated random field value: " << wArrays_[0][2] << std::endl;
         
         // If flexible, do mixing of histories for unit cell parameters
         if (isFlexible_) {
            for (int m = 0; m < nParameter ; ++m) {
               wCpArrays_[m] = cellParamHists_[0][m];
               dCpArrays_[m] = stressHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < nParameter ; ++m) {
                  wCpArrays_[m] += coeffs_[i]*( cellParamHists_[i+1][m] -
                                                cellParamHists_[0][m]);
                  dCpArrays_[m] += coeffs_[i]*( stressHists_[i+1][m] -
                                                stressHists_[0][m]);
               }
            }
            parameters_.clear();
            for (int m = 0; m < nParameter ; ++m) {
               parameters_.append(wCpArrays_[m] + lambda_ * dCpArrays_[m]);
            }
            system().setUnitCell(parameters_);
         }
      }
   }

   template <int D>
   bool AmIterator<D>::isCanonical()
   {
      // check ensemble of all polymers
      for (int i = 0; i < system().mixture().nPolymer(); ++i) {
         if (system().mixture().polymer(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // check ensemble of all solvents
      for (int i = 0; i < system().mixture().nSolvent(); ++i) {
         if (system().mixture().solvent(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // returns true if false was never returned
      return true;

   }

   template <int D>
   void AmIterator<D>::cleanUp()
   {
      // Deallocate arrays used in coefficient calculation.
      if (U_.isAllocated()) {
         U_.deallocate();
      }
      if (coeffs_.isAllocated()) {
         coeffs_.deallocate();
      }
      if (v_.isAllocated()) {
         v_.deallocate();
      }

   }

}
}
#endif
