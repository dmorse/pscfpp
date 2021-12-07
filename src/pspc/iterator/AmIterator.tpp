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
      devHists_.allocate(maxHist_+1);
      omHists_.allocate(maxHist_+1);

      if (isFlexible_) {
         devCpHists_.allocate(maxHist_+1);
         CpHists_.allocate(maxHist_+1);
      }

      int nMonomer = system().mixture().nMonomer();
      wArrays_.allocate(nMonomer);
      dArrays_.allocate(nMonomer);
      tempDev.allocate(nMonomer);

      int nBasis = system().basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         wArrays_[i].allocate(nBasis - 1);
         dArrays_[i].allocate(nBasis - 1);
         tempDev[i].allocate(nBasis - 1);
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

      FieldIo<D>& fieldIo = system().fieldIo();

      // Solve MDE for initial state
      solverTimer.start();
      system().mixture().compute(system().wFieldsRGrid(),
                                    system().cFieldsRGrid());
      now = Timer::now();
      solverTimer.stop(now);

      // Convert c fields from RGrid to Basis
      convertTimer.start(now);
      fieldIo.convertRGridToBasis(system().cFieldsRGrid(),
                                  system().cFields());
      now = Timer::now();
      convertTimer.stop(now);

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
         computeDeviation();

         // Test for convergence
         done = isConverged();

         if (done) {

            updateTimer.stop();
            Log::file() << "----------CONVERGED----------"<< std::endl;

            // Output timing results
            double updateTime = updateTimer.time();
            double convertTime = convertTimer.time();
            double solverTime = solverTimer.time();
            double stressTime = 0.0;
            double totalTime = updateTime + convertTime + solverTime;
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
               for (int m=0; m < system().unitCell().nParameter(); ++m){
                  Log::file() << "Stress  "<< m << "   = "
                              << Dbl(system().mixture().stress(m)) 
                              << "\n";
               }
               Log::file() << "\n";
            }

            // Successful completion (i.e., converged within tolerance)
            cleanUp();
            std::cout << "\n Iterations = " << itr << std::endl;
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
               if (nHist_ > 0) {
                  invertMatrix_.deallocate();
                  coeffs_.deallocate();
                  vM_.deallocate();
               }
            }
            now = Timer::now();
            updateTimer.stop(now);

            // Convert wFields from Basis to RGrid
            convertTimer.start(now);
            fieldIo.convertBasisToRGrid(system().wFields(),
                                        system().wFieldsRGrid());
            now = Timer::now();
            convertTimer.stop(now);

            // Solve MDE
            solverTimer.start(now);
            system().mixture().compute(system().wFieldsRGrid(),
                                       system().cFieldsRGrid());
            now = Timer::now();
            solverTimer.stop(now);

            // Compute stress if needed
            if (isFlexible_){
               stressTimer.start(now);
               system().mixture().computeStress();
               now = Timer::now();
               stressTimer.stop(now);
            }

            // Transform computed cFields from RGrid to Basis
            convertTimer.start(now);
            fieldIo.convertRGridToBasis(system().cFieldsRGrid(),
                                        system().cFields());
            now = Timer::now();
            convertTimer.stop(now);
         }

      }

      // Failure: iteration counter itr reached maxItr without converging
      cleanUp();
      return 1;
   }

   template <int D>
   void AmIterator<D>::computeDeviation()
   {

      omHists_.append(system().wFields());

      if (isFlexible_)
         CpHists_.append(system().unitCell().parameters());

      for (int i = 0 ; i < system().mixture().nMonomer(); ++i) {
         for (int j = 0; j < system().basis().nBasis() - 1; ++j) {
            tempDev[i][j] = 0;
         }
      }

      DArray<double> temp;
      temp.allocate(system().basis().nBasis() - 1);

      #if 0
      for (int i = 0; i < system().mixture().nMonomer(); ++i) {

         for (int j = 0; j < system().basis().nBasis() - 1; ++j) {
            temp[j] = 0;
         }

         for (int j = 0; j < system().mixture().nMonomer(); ++j) {
            for (int k = 0; k < system().basis().nBasis() - 1; ++k) {
               tempDev[i][k] += system().interaction().chi(i,j) *
                              system().cField(j)[k + 1];
               temp[k] += system().wField(j)[k + 1];
            }
         }

         for (int k = 0; k < system().basis().nBasis() - 1; ++k) {
            tempDev[i][k] += ((temp[k] / system().mixture().nMonomer())
                             - system().wField(i)[k + 1]);
         }
      }
      #endif

      for (int i = 0; i < system().mixture().nMonomer(); ++i) {
         for (int j = 0; j < system().mixture().nMonomer(); ++j) {
            for (int k = 0; k < system().basis().nBasis() - 1; ++k) {
               tempDev[i][k] +=( (system().interaction().chi(i,j)*system().cField(j)[k + 1])
                               - (system().interaction().idemp(i,j)*system().wField(j)[k + 1]) );
            }
         }
      }

      devHists_.append(tempDev);

      if (isFlexible_){
         FArray<double, 6 > tempCp;
         for (int i = 0; i < system().unitCell().nParameter() ; i++){
            tempCp [i] = -((system().mixture()).stress(i));
         }
         devCpHists_.append(tempCp);
      }
   }

   template <int D>
   bool AmIterator<D>::isConverged()
   {
      double error;

      #if 0
      // Error as defined in Matsen's Papers
      double dError = 0;
      double wError = 0;
      for ( int i = 0; i < system().mixture().nMonomer(); i++) {
         for ( int j = 0; j < system().basis().nBasis() - 1; j++) {
            dError += devHists_[0][i][j] * devHists_[0][i][j];

            //the extra shift is due to the zero indice coefficient being
            //exactly known
            wError += system().wField(i)[j+1] * system().wField(i)[j+1];
         }
      }

      if (isFlexible_){
         for ( int i = 0; i < system().unitCell().nParameter() ; i++) {
            dError += devCpHists_[0][i] *  devCpHists_[0][i];
            wError += system().unitCell().parameters()[i]
                     *system().unitCell().parameters()[i];
         }
      }
      Log::file() << " dError :" << Dbl(dError)<<std::endl;
      Log::file() << " wError :" << Dbl(wError)<<std::endl;
      error = sqrt(dError / wError);
      #endif

      // Error by Max Residuals
      double temp1 = 0;
      double temp2 = 0;
      for ( int i = 0; i < system().mixture().nMonomer(); i++) {
         for ( int j = 0; j < system().basis().nBasis() - 1; j++) {
            if (temp1 < fabs (devHists_[0][i][j]))
                temp1 = fabs (devHists_[0][i][j]);
         }
      }
      Log::file() << "SCF Error   = " << Dbl(temp1) << std::endl;
      error = temp1;

      if (isFlexible_){
         for ( int i = 0; i < system().unitCell().nParameter() ; i++) {
            if (temp2 < fabs (devCpHists_[0][i])) {
                temp2 = fabs (devCpHists_[0][i]);
            }
         }
         // Output current stress values
         for (int m=0;  m < system().unitCell().nParameter() ; ++m){
            Log::file() << "Stress  "<< m << "   = "
                        << Dbl(system().mixture().stress(m)) <<"\n";
         }
         error = (temp1>(100*temp2)) ? temp1 : (100*temp2);
         // 100 is chose as stress rescale factor
         // TODO: Separate SCF and stress tolerance limits
      }
      Log::file() << "Error       = " << Dbl(error) << std::endl;

      // Output current unit cell parameter values
      if (isFlexible_){
         for (int m=0; m < system().unitCell().nParameter() ; ++m){
               Log::file() << "Parameter " << m << " = "
                           << Dbl(system().unitCell().parameters()[m])
                           << "\n";
         }
      }

      // Check if total error is below tolerance
      if (error < epsilon_) {
         return true;
      } else {
         return false;
      }
   }

   template <int D>
   void AmIterator<D>::minimizeCoeff(int itr)
   {
      if (itr == 1) {
         //do nothing
      } else {

         int nMonomer = system().mixture().nMonomer();
         int nParameter = system().unitCell().nParameter();
         int nBasis = system().basis().nBasis();
         double elm, elm_cp;

         for (int i = 0; i < nHist_; ++i) {
            for (int j = i; j < nHist_; ++j) {

               invertMatrix_(i,j) = 0;
               for (int k = 0; k < nMonomer; ++k) {
                  elm = 0;
                  for (int l = 0; l < nBasis - 1; ++l) {
                     elm +=
                            ((devHists_[0][k][l] - devHists_[i+1][k][l])*
                             (devHists_[0][k][l] - devHists_[j+1][k][l]));
                  }
                  invertMatrix_(i,j) += elm;
               }

               if (isFlexible_){
                  elm_cp = 0;
                  for (int m = 0; m < nParameter ; ++m){
                     elm_cp += ((devCpHists_[0][m] - devCpHists_[i+1][m])*
                                (devCpHists_[0][m] - devCpHists_[j+1][m]));
                  }
                  invertMatrix_(i,j) += elm_cp;
               }
               invertMatrix_(j,i) = invertMatrix_(i,j);
            }

            vM_[i] = 0;
            for (int j = 0; j < nMonomer; ++j) {
               for (int k = 0; k < nBasis - 1; ++k) {
                  vM_[i] += ( (devHists_[0][j][k] - devHists_[i+1][j][k]) *
                               devHists_[0][j][k] );
               }
            }

            if (isFlexible_){
               elm_cp = 0;
               for (int m = 0; m < nParameter ; ++m){
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
      UnitCell<D>& unitCell = system().unitCell();
      Mixture<D>&  mixture = system().mixture();

      if (itr == 1) {
         for (int i = 0; i < mixture.nMonomer(); ++i) {
            for (int j = 0; j < system().basis().nBasis() - 1; ++j) {
               system().wField(i)[j+1]
                      = omHists_[0][i][j+1] + lambda_*devHists_[0][i][j];
            }
         }

         if (isFlexible_){
            parameters_.clear();
            for (int m = 0; m < unitCell.nParameter() ; ++m){
               parameters_.append(CpHists_[0][m]
                              + lambda_* devCpHists_[0][m]);

            }
            unitCell.setParameters(parameters_);
            mixture.setupUnitCell(unitCell);
            system().basis().update();
         }

      } else {
         for (int j = 0; j < mixture.nMonomer(); ++j) {
            for (int k = 0; k < system().basis().nBasis() - 1; ++k) {
               wArrays_[j][k] = omHists_[0][j][k + 1];
               dArrays_[j][k] = devHists_[0][j][k];
            }
         }
         for (int i = 0; i < nHist_; ++i) {
            for (int j = 0; j < mixture.nMonomer(); ++j) {
               for (int k = 0; k < system().basis().nBasis() - 1; ++k) {
                  wArrays_[j][k] += coeffs_[i] * ( omHists_[i+1][j][k+1] -
                                                   omHists_[0][j][k+1] );
                  dArrays_[j][k] += coeffs_[i] * ( devHists_[i+1][j][k] -
                                                   devHists_[0][j][k] );
               }
            }
         }
         for (int i = 0; i < mixture.nMonomer(); ++i) {
            for (int j = 0; j < system().basis().nBasis() - 1; ++j) {
               system().wField(i)[j+1] = wArrays_[i][j]
                                         + lambda_ * dArrays_[i][j];
            }
         }
         if (isFlexible_){
            for (int m = 0; m < unitCell.nParameter() ; ++m){
               wCpArrays_[m] = CpHists_[0][m];
               dCpArrays_[m] = devCpHists_[0][m];
            }
            for (int i = 0; i < nHist_; ++i) {
               for (int m = 0; m < unitCell.nParameter() ; ++m) {
                  wCpArrays_[m] += coeffs_[i]*( CpHists_[i+1][m] -
                                                CpHists_[0][m]);
                  dCpArrays_[m] += coeffs_[i]*( devCpHists_[i+1][m] -
                                                devCpHists_[0][m]);
               }
            }
            parameters_.clear();
            for (int m = 0; m < unitCell.nParameter() ; ++m){
               parameters_.append(wCpArrays_[m] + lambda_ * dCpArrays_[m]);
            }
            unitCell.setParameters(parameters_);
            mixture.setupUnitCell(unitCell);
            system().basis().update();
         }
      }
   }

   template <int D>
   void AmIterator<D>::cleanUp()
   {
      // Deallocate allocated arrays and matrices, if allocated.
      if (invertMatrix_.isAllocated()) {
         invertMatrix_.deallocate();
      }
      if (coeffs_.isAllocated()) {
         coeffs_.deallocate();
      }
      if (vM_.isAllocated()) {
         vM_.deallocate();
      }

      
   }

}
}
#endif
