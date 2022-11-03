#ifndef PSPC_MC_SIMULATOR_TPP
#define PSPC_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"
#include <pspc/System.h>
#include <pspc/mcmove/McMoveFactory.h>

#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>

#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
   : Manager< McMove<D> >(),
     systemPtr_(&system),
     randomPtr_(&system.random())
   {  setClassName("McSimulator"); }

   /*
   * Destructor.
   */
   template <int D>
   McSimulator<D>::~McSimulator()
   {}

   /*
   * Return a pointer to a new McMoveFactory object.
   */
   template <int D>
   Factory< McMove<D> >* McSimulator<D>::newDefaultFactory() const
   {  return new McMoveFactory<D>(*systemPtr_); }

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McSimulator<D>::readParameters(std::istream &in)
   {
      Manager< McMove<D> >::readParameters(in);

      if (size() > 0) {

         // Allocate and store probabilities
         probabilities_.allocate(size());
         double  totalProbability = 0.0;
         int     iMove;
         for (iMove = 0; iMove < size(); ++iMove) {
            probabilities_[iMove] = (*this)[iMove].probability();
            totalProbability += probabilities_[iMove];
         }
   
         // Allocate and store and normalize probabilities
         for (iMove = 0; iMove < size(); ++iMove) {
            probabilities_[iMove] = probabilities_[iMove]/totalProbability;
            (*this)[iMove].setProbability(probabilities_[iMove]);
         }

      }

      // Allocate projected chi matrix chiP_ and associated arrays
      const int nMonomer = system().mixture().nMonomer();
      chiP_.allocate(nMonomer, nMonomer);
      chiEvals_.allocate(nMonomer);
      chiEvecs_.allocate(nMonomer, nMonomer);

      analyzeChi();
   }

   /*
   * Initialize just prior to a run.
   */
   template <int D>
   void McSimulator<D>::setup()
   {

      // Allocate mcState_, if necessary.
      if (!mcState_.isAllocated) {
         const int nMonomer = system().mixture().nMonomer();
         const IntVec<D> dimensions = system().domain().mesh().dimensions();
         mcState_.allocate(nMonomer, dimensions);
      }

      // Call setup fucntion of each McMove
      for (int iMove = 0; iMove < size(); ++iMove) {
         (*this)[iMove].setup();
      }

   }

   /*
   * Choose a McMove at random.
   */
   template <int D>
   McMove<D>& McSimulator<D>::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   template <int D>
   void McSimulator<D>::output()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void McSimulator<D>::simulate(int nStep)
   {

      setup();
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for (int iStep = 0; iStep < nStep; ++iStep) {

         // Choose and attempt an McMove
         chooseMove().move();

      }
      timer.stop();
      double time = timer.time();

      // Output results of move statistics to files
      output();

      // Output time for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep         " << nStep << std::endl;
      Log::file() << "run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep  " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << std::endl;

      // Print McMove acceptance statistics
      long attempt;
      long accept;
      using namespace std;
      Log::file() << "Move Statistics:" << endl << endl;
      Log::file() << setw(32) << left <<  "Move Name"
           << setw(12) << right << "Attempted"
           << setw(12) << right << "Accepted"
           << setw(15) << right << "AcceptRate"
           << endl;
      int nMove = size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = (*this)[iMove].nAttempt();
         accept  = (*this)[iMove].nAccept();
         Log::file() << setw(32) << left
              << (*this)[iMove].className()
              << setw(12) << right << attempt
              << setw(12) << accept
              << setw(15) << fixed << setprecision(6)
              << ( attempt == 0 ? 0.0 : double(accept)/double(attempt) )
              << endl;
      }
      Log::file() << endl;

   }

   /*
   * Save the current Monte-Carlo state.
   *
   * Used before attempting a Monte-Carlo move.
   */
   template <int D>
   void McSimulator<D>::saveMcState()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(!mcState_.hasData);

      int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         mcState_.w[i] = system().w().rgrid(i);
      }
      mcState_.mcHamiltonian  = mcHamiltonian_;
      mcState_.hasData = true;
   }

   /*
   * Restore a saved Monte-Carlo state.
   *
   * Used when an attempted Monte-Carlo move is rejected.
   */
   template <int D>
   void McSimulator<D>::restoreMcState()
   {
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(mcState_.hasData);

      system().setWRGrid(mcState_.w); 
      mcHamiltonian_ = mcState_.mcHamiltonian;
      mcState_.hasData = false;
      hasMcHamiltonian_ = true;
   }

   /*
   * Compute Monte Carlo Hamiltonian.
   */
   template <int D>
   void McSimulator<D>::computeMcHamiltonian()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().hasCFields());

      // Compute free energy if necessary
      if (!system().hasFreeEnergy()) {
         system().computeFreeEnergy();
      }

      //double h = system().fHelmholtz();

      hasMcHamiltonian_ = true;
   }

   template <int D>
   void McSimulator<D>::analyzeChi()
   {
      const int nMonomer = system().mixture().nMonomer();
      DMatrix<double> const & chi = system().interaction().chi();
      double d = 1.0/double(nMonomer);
      int i, j, k;

      // Compute projection matrix P
      DMatrix<double> P;
      P.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            P(i,j) = -d;
         }
         P(i,i) += 1.0;
      }

      // Compute T = chi*P (temporary matrix)
      DMatrix<double> T;
      T.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            T(i, j) = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               T(i,j) += chi(i,k)*P(k,j);
            }
         }
      }

      // Allocate GSL matrix to hold elements of chiP
      gsl_matrix* A = gsl_matrix_alloc(nMonomer, nMonomer);

      // Compute chiP = = P*chi*P = P*T
      DMatrix<double> chiP;
      chiP.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            chiP(i, j) = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               chiP(i,j) += P(i,k)*T(k,j);
            }
            gsl_matrix_set(A, i, j, chiP(i, j));
         }
      }

      // Compute eigenvalues and eigenvectors of chiP
      gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(nMonomer);
      gsl_vector* Avals = gsl_vector_alloc(nMonomer);
      gsl_matrix* Avecs = gsl_matrix_alloc(nMonomer, nMonomer);
      int error;
      error = gsl_eigen_symmv(A, Avals, Avecs, work);

      // Copy eigenpairs with non-null eigenvalues
      int nNull =  0;
      int iNull = -1;
      k = 0;
      double val;
      for (i = 0; i < nMonomer; ++i) {
         val = gsl_vector_get(Avals, i);
         if (abs(val) < 1.0E-8) {
            ++nNull;
            iNull = i;
         } else {
            chiEvals_[k] = val;
            for (j = 0; j < nMonomer; ++j) {
               chiEvecs_(k, j) = gsl_matrix_get(Avecs, j, i);
            }
            if (chiEvecs_(k, 0) < 0.0) {
               for (j = 0; j < nMonomer; ++j) {
                  chiEvecs_(k, j) = -chiEvecs_(k, j);
               }
            }
            ++k;
         }
      }
      UTIL_CHECK(nNull == 1);
      UTIL_CHECK(iNull >= 0);
    
      // Copy eigenpair with null eigenvalue 
      chiEvals_[nMonomer-1] = 0.0;
      for (i = 0; i < nMonomer; ++i) {
         chiEvecs_(nMonomer-1, i) = gsl_matrix_get(Avecs, i, iNull);
      }

      #if 0
      // Debugging 
      for (i = 0; i < nMonomer; ++i) {
         Log::file() << "Eigenpair " << i << "\n";
         Log::file() << "value  =  " << chiEvals_[i] << "\n";
         Log::file() << "vector = [ "; 
         for (j = 0; j < nMonomer; ++j) {
            Log::file() << chiEvecs_(i, j) << "   ";
         }
         Log::file() << "]\n";
      }
      #endif

   }
}
}
#endif
