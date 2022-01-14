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
      scaleStress_ = 10.0; // default stress scaling
      errorType_ = "normResid"; // default type of error
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "errorType", errorType_);
      readOptional(in, "scaleStress", scaleStress_);

      if (!(errorType_ == "normResid" || errorType_ == "maxResid")) {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }
      
   }

   /*
   * Setup and allocate memory required by iterator.
   */
   template <int D>
   void AmIterator<D>::setup()
   {
      // Determine length of residual basis vectors
      nResid_ = system().mixture().nMonomer();
      if (isFlexible_) {
         nResid_ += system().unitCell().nParameter(); 
      }

      // Determine how to treat homogeneous basis function coefficients
      if (isCanonical()) {
         // Directly calculate them
         shift_ = 1;
         isCanonical_ = true;
      } else {
         // Iteratively calculate them
         shift_ = 0;
         isCanonical_ = false;
      }

      // Allocate ring buffers
      resHists_.allocate(maxHist_+1);
      wHists_.allocate(maxHist_+1);

      if (isFlexible_) {
         stressHists_.allocate(maxHist_+1);
         cellParamHists_.allocate(maxHist_+1);
      }

      // Allocate outer arrays used in iteration
      int nMonomer = system().mixture().nMonomer();
      wArrays_.allocate(nMonomer);
      dArrays_.allocate(nMonomer);
      resArrays_.allocate(nResid_);
      
      // Allocate inner arrays with number of basis functions
      int nBasis = system().basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         wArrays_[i].allocate(nBasis);
         dArrays_[i].allocate(nBasis);
      }

      // Allocate inner residual arrays
      for (int i = 0; i < nResid_; ++i) {
         resArrays_[i].allocate(nElem(i));
      }

      // Allocate arrays/matrices used in coefficient calculation
      U_.allocate(maxHist_,maxHist_);
      v_.allocate(maxHist_);
      coeffs_.allocate(maxHist_);
   }

   /*
   * Solve iteratively.
   */
   template <int D>
   int AmIterator<D>::solve()
   {

      // Preconditions:
      UTIL_CHECK(system().hasWFields());

      // Timers for timing iteration components
      // An idea: timer "dictionary" to make iterating over timers easy
      // Automating output with flexible numbers of timers
      // DArray< Timer > timers;
      // DArray< std::string > timernames;

      Timer timerMDE;
      Timer timerStress;
      Timer timerAM;
      Timer timerResid;
      Timer timerConverged;
      Timer timerCoeff;
      Timer timerOmega;
      Timer timerTotal;
      

      // Start overall timer 
      timerTotal.start();
      
      // Solve MDE for initial state
      timerMDE.start();
      system().compute();
      timerMDE.stop();

      // Compute initial stress if needed
      if (isFlexible_) {
         timerStress.start();
         system().mixture().computeStress();
         timerStress.stop();
      }

      // Iterative loop
      bool done;
      for (int itr = 0; itr < maxItr_; ++itr) {

         timerAM.start();

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<" Iteration  "<<itr<<std::endl;

         if (itr < maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr+1);
            nHist_ = itr;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }
         timerResid.start();
         computeResidual();
         timerResid.stop();

         // Test for convergence
         timerConverged.start();
         done = isConverged();
         timerConverged.stop();

         if (done) {
            // Stop timers
            timerAM.stop();
            timerTotal.stop();

            Log::file() << "----------CONVERGED----------"<< std::endl;

            // Output timing results
            Log::file() << "\n";
            Log::file() << "Iterator times contributions:\n";
            Log::file() << "\n";
            Log::file() << "MDE solution:         " 
                        << timerMDE.time()  << " s,  "
                        << timerMDE.time()/timerTotal.time() << "\n";
            if (isFlexible_) {
               Log::file() << "stress computation:   " 
                           << timerStress.time()  << " s,  "
                           << timerStress.time()/timerTotal.time() << "\n";
            }
            Log::file() << "residual computation: "  
                        << timerResid.time()  << " s,  "
                        << timerResid.time()/timerTotal.time() << "\n";
            Log::file() << "mixing coefficients:  "  
                        << timerCoeff.time()  << " s,  "
                        << timerCoeff.time()/timerTotal.time() << "\n";
            Log::file() << "checking convergence: "  
                        << timerConverged.time()  << " s,  "
                        << timerConverged.time()/timerTotal.time() << "\n";
            Log::file() << "updating field guess: "  
                        << timerOmega.time()  << " s,  "
                        << timerOmega.time()/timerTotal.time() << "\n";
            Log::file() << "total time:           "  
                        << timerTotal.time()   << " s  ";
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
            // Determine optimal linear combination coefficients for 
            // building the updated field guess
            timerCoeff.start();
            minimizeCoeff();
            timerCoeff.stop();
            // Build the updated field
            timerOmega.start();
            buildOmega();
            timerOmega.stop();

            timerAM.stop();

            // Solve MDE
            timerMDE.start();
            system().compute();
            timerMDE.stop();

            // Compute stress if needed
            if (isFlexible_) {
               timerStress.start();
               system().mixture().computeStress();
               timerStress.stop();
            }

         }

      }
      // Failure: iteration counter itr reached maxItr without converging
      timerTotal.stop();
      
      Log::file() << "Iterator failed to converge before the maximum number of iterations.\n";
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
      for (int i = 0 ; i < nResid_; ++i) {
         for (int k = 0; k < nElem(i); ++k) {
            resArrays_[i][k] = 0.0;
         }
      }
      
      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = shift_; k < nBasis; ++k) {
               resArrays_[i][k] +=
                  system().interaction().chi(i,j)*system().cField(j)[k] -
                  system().interaction().idemp(i,j)*system().wField(j)[k];
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (!isCanonical_) {
         for (int i = 0; i < nMonomer; ++i) {
            resArrays_[i][0] -= 1.0/system().interaction().sum_inv();
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible_) {
         FArray<double, 6 > tempCp;
         for (int i = 0; i < nParameter ; i++) {
            tempCp[i] = -((system().mixture()).stress(i));
         }
         stressHists_.append(tempCp);

         for (int m=0, i = nMonomer; m < nParameter ; ++m, ++i) {
            resArrays_[i][0] = scaleStress_ * fabs( stressHists_[0][m] );
         }
      }
      
      // Store residuals in residual history ringbuffer
      resHists_.append(resArrays_);

      return;
   }

   template <int D>
   bool AmIterator<D>::isConverged()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();

      // Find max residual
      double maxSCF = 0.0, maxStress = 0.0, maxRes = 0.0;

      for (int i = 0; i < nResid_; ++i) {
         for (int k = 0; k < nElem(i); ++k) {
            double currRes = fabs(resHists_[0][i][k]);
            if (i < nMonomer && maxSCF < currRes) {
               maxSCF = currRes;
            } else
            if (i >= nMonomer && maxStress < currRes) {
               maxStress = currRes;
            }
         }
      }

      Log::file() << "SCF Error   = " << Dbl(maxSCF) << std::endl;
      maxRes = maxSCF;

      // Error by max Stress residual
      if (isFlexible_) { 
         // Check if stress residual is greater than SCF
         if (maxStress > maxRes) {
            maxRes = maxStress;
         }
         // Output stress values
         for (int m=0;  m < nParameter ; ++m) {
            Log::file() << "Stress  "<< m << "   = "
                        << Dbl(system().mixture().stress(m)) <<"\n";
         }
      }
      Log::file() << "Max Residual = " << Dbl(maxRes) << std::endl;

      // Output current unit cell parameter values
      if (isFlexible_) {
         for (int m=0; m < nParameter ; ++m) {
               Log::file() << "Parameter " << m << " = "
                           << Dbl(system().unitCell().parameters()[m])
                           << "\n";
         }
      }

      // Error by norm of residual vector
      double normResSq = 0.0, normRes = 0.0;
      for (int i = 0; i < nResid_; ++i) {
         for (int k = 0; k < nElem(i); ++k) {
            normResSq += resHists_[0][i][k] * resHists_[0][i][k];
         }
      }
      normRes = sqrt(normResSq);

      Log::file() << "Residual Norm = " << Dbl(normRes) << std::endl;

      // Check if total error is below tolerance
      if (errorType_ == "normResid") {
         return normRes < epsilon_;
      } else 
      if (errorType_ == "maxResid") {
         return maxRes < epsilon_;
      } else {
         UTIL_THROW("Invalid iterator error type in parameter file.");
         return false;
      }

   }

   template <int D>
   void AmIterator<D>::minimizeCoeff()
   {
      // Initialize matrix and vector of residual dot products
      // if this is the first iteration
      if (nHist_ == 0) {
         for (int m = 0; m < maxHist_; ++m) {
            v_[m] = 0;
            coeffs_[m] = 0;
            for (int n = 0; n < maxHist_; ++n) {
               U_(m,n) = 0;
            }
         }
         return;
      }

      // Update matrix U by shifting elements diagonally
      for (int m = maxHist_-1; m > 0; --m) {
         for (int n = maxHist_-1; n > 0; --n) {
            U_(m,n) = U_(m-1,n-1); 
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist_; ++m) {
         // Compute basis vector dot product
         double dotprod = 0;
         for (int i = 0; i < nResid_; ++i) {
            for (int k = 0; k < nElem(i); ++k) {
               dotprod += (resHists_[m][i][k] - resHists_[m+1][i][k]) *
                          (resHists_[0][i][k] - resHists_[1][i][k]);
            }
         }
         U_(m,0) = dotprod;
         U_(0,m) = dotprod;
      }

      // Compute v vector, as described in Arora 2017.
      for (int m = 0; m < nHist_; ++m) {
         v_[m] = 0;
         // Dot product of residual vectors
         for (int i = 0; i < nResid_; ++i) {
            for (int k = 0; k < nElem(i); ++k) {
               v_[m] += resHists_[0][i][k] * 
                        (resHists_[m][i][k] - resHists_[m+1][i][k]);
            }
         }
      }
      
      // Solve matrix equation problem to get coefficients to minimize
      // the norm of the residual vector
      if (nHist_ == 1) {
         // Solve explicitly for coefficient
         coeffs_[0] = v_[0] / U_(0,0);
      } else
      if (nHist_ < maxHist_) { 
         // Create temporary smaller version of U_, v_, coeffs_
         // this is done to avoid reallocating U_ with each iteration.
         DMatrix<double> tempU;
         DArray<double> tempv,tempcoeffs;
         tempU.allocate(nHist_,nHist_);
         tempv.allocate(nHist_);
         tempcoeffs.allocate(nHist_);
         for (int i = 0; i < nHist_; ++i) {
            tempv[i] = v_[i];
            for (int j = 0; j < nHist_; ++j) {
               tempU(i,j) = U_(i,j);
            }
         }
         // Solve matrix equation
         LuSolver solver;
         solver.allocate(nHist_);
         solver.computeLU(tempU);
         solver.solve(tempv,tempcoeffs);
         // Transfer solution to full-sized member variable
         for (int i = 0; i < nHist_; ++i) {
            coeffs_[i] = tempcoeffs[i];
         }
      } else 
      if (nHist_ == maxHist_) {
         LuSolver solver;
         solver.allocate(maxHist_);
         solver.computeLU(U_);
         solver.solve(v_, coeffs_);
      }

      return;
   }

   template <int D>
   void AmIterator<D>::buildOmega()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nParameter = system().unitCell().nParameter();
      const int nBasis = system().basis().nBasis();

      if (nHist_ == 0) { // if 0 historical solutions

         // Update omega field with SCF residuals
         for (int i = 0; i < nMonomer; ++i) {
            for (int k = shift_; k < nBasis; ++k) {
               wArrays_[i][k]
                      = wHists_[0][i][k] + lambda_*resHists_[0][i][k];
            }
         }
         // If canonical, set the homogeneous components explicitly
         if (isCanonical_) {
            for (int i = 0; i < nMonomer; ++i) {
               dArrays_[i][0] = 0.0;
               wArrays_[i][0] = 0.0;
               for (int j = 0; j < nMonomer; ++j) {
                  wArrays_[i][0] += 
                     system().interaction().chi(i,j) * system().cField(j)[0];
               }
            }
         }
         system().setWBasis(wArrays_);

         // Update unit cell parameters with stress residuals
         if (isFlexible_) {
            parameters_.clear();
            for (int m = 0; m < nParameter ; ++m) {
               parameters_.append(cellParamHists_[0][m]
                              + lambda_ * stressHists_[0][m]);

            }
            system().setUnitCell(parameters_);
         }

      } else { // if at least 1 historical result
         // Contribution of the last solution
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = shift_; k < nBasis; ++k) {
               wArrays_[j][k] = wHists_[0][j][k];
               dArrays_[j][k] = resHists_[0][j][k];
            }
         }
         // Mixing in historical solutions
         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < nMonomer; ++j) {
               for (int k = shift_; k < nBasis; ++k) {
                  wArrays_[j][k] += coeffs_[i] * ( wHists_[i+1][j][k] -
                                                   wHists_[i][j][k] );
                  dArrays_[j][k] += coeffs_[i] * ( resHists_[i+1][j][k] -
                                                   resHists_[i][j][k] );
               }
            }
         }
         // If canonical ensemble, set the homogeneous components explicitly
         if (isCanonical_) {
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
         
         // If flexible, do mixing of histories for unit cell parameters
         if (isFlexible_) {
            for (int m = 0; m < nParameter ; ++m) {
               wCpArrays_[m] = cellParamHists_[0][m];
               dCpArrays_[m] = stressHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < nParameter ; ++m) {
                  wCpArrays_[m] += coeffs_[i] * ( cellParamHists_[i+1][m] -
                                                  cellParamHists_[i][m] );
                  dCpArrays_[m] += coeffs_[i] * ( stressHists_[i+1][m] -
                                                  stressHists_[i][m] );
               }
            }
            parameters_.clear();
            for (int m = 0; m < nParameter ; ++m) {
               parameters_.append(wCpArrays_[m] + lambda_ * dCpArrays_[m]);
            }
            system().setUnitCell(parameters_);
         }
      }
      return;
   }

   template <int D>
   bool AmIterator<D>::isCanonical()
   {
      // Check ensemble of all polymers
      for (int i = 0; i < system().mixture().nPolymer(); ++i) {
         if (system().mixture().polymer(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Check ensemble of all solvents
      for (int i = 0; i < system().mixture().nSolvent(); ++i) {
         if (system().mixture().solvent(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Returns true if false was never returned
      return true;
   }

   template <int D>
   int AmIterator<D>::nElem(int i)
   {
      const int nMonomer = system().mixture().nMonomer();
      int nBasis = system().basis().nBasis();

      if (i < nMonomer) {
         return nBasis;
      } else {
         return 1;
      }
   }

   template <int D>
   void AmIterator<D>::cleanUp()
   {
      // Clear ring buffers
      resHists_.clear();
      wHists_.clear();
      stressHists_.clear();
      cellParamHists_.clear();

      return;
   }

}
}
#endif
