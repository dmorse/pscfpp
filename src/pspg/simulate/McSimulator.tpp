#ifndef PSPG_MC_SIMULATOR_TPP
#define PSPG_MC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSimulator.h"
#include <pspg/System.h>
#include <pspg/simulate/mcmove/McMoveFactory.h>
#include <pspg/simulate/analyzer/AnalyzerFactory.h>
#include <pspg/simulate/trajectory/TrajectoryReader.h>
#include <pspg/simulate/trajectory/TrajectoryReaderFactory.h>
#include <pspg/compressor/Compressor.h>
#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>
#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McSimulator<D>::McSimulator(System<D>& system)
   : mcMoveManager_(*this, system),
     analyzerManager_(*this, system),
     trajectoryReaderFactoryPtr_(0),
     random_(),
     systemPtr_(&system),
     hasMcHamiltonian_(false),
     hasWC_(false)
   { 
      setClassName("McSimulator");   
      trajectoryReaderFactoryPtr_ = new TrajectoryReaderFactory<D>(system);
   }

   /*
   * Destructor.
   */
   template <int D>
   McSimulator<D>::~McSimulator()
   {}

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McSimulator<D>::readParameters(std::istream &in)
   {
      // Read block of mc move data inside
      readParamComposite(in, mcMoveManager_);
      // Read block of analyzer
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);
      // Allocate projected chi matrix chiP_ and associated arrays
      const int nMonomer = system().mixture().nMonomer();
      chiP_.allocate(nMonomer, nMonomer);
      chiEvals_.allocate(nMonomer);
      chiEvecs_.allocate(nMonomer, nMonomer);
      wc_.allocate(nMonomer);
      const int meshSize = system().domain().mesh().size();
      for (int i = 0; i < nMonomer; ++i) {
            wc_[i].allocate(meshSize);
      }
      analyzeChi();
   }

   /*
   * Initialize just prior to a run.
   */
   template <int D>
   void McSimulator<D>::setup()
   {  
      UTIL_CHECK(system().w().hasData());
      // Allocate mcState_, if necessary.
      if (!mcState_.isAllocated) {
         const int nMonomer = system().mixture().nMonomer();
         int meshSize = system().domain().mesh().size();
         mcState_.allocate(nMonomer, meshSize);
      }
   
      // Eigenanalysis of the projected chi matrix.
      analyzeChi();
      // Compute field components and MC Hamiltonian for initial state
      system().compute();
      computeWC();
      computeMcHamiltonian();
      mcMoveManager_.setup();
      if (analyzerManager_.size() > 0){
         analyzerManager_.setup();
      }
   }

   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void McSimulator<D>::simulate(int nStep)
   {
      UTIL_CHECK(mcMoveManager_.size() > 0);

      setup();
      Log::file() << std::endl;

      // Main Monte Carlo loop
      Timer timer;
      timer.start();
      for (int iStep = 0; iStep < nStep; ++iStep) {
         // Analysis (if any)
         if (Analyzer<D>::baseInterval != 0) {
            if (iStep % Analyzer<D>::baseInterval == 0) {
               if (analyzerManager_.size() > 0) {
                  iStep_ = iStep;
                  analyzerManager_.sample(iStep);
               }
            }
         }
         // Choose and attempt an McMove
         mcMoveManager_.chooseMove().move();
      }
      timer.stop();
      double time = timer.time();

      // Output results of move statistics to files
      mcMoveManager_.output();
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }
      
      // Output how many times MDE has been solved for the simulation run
      Log::file() << std::endl;
      Log::file() << "MDE counter   " << system().compressor().counterMDE()<< std::endl;
      Log::file() << std::endl;
      
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
      int nMove = mcMoveManager_.size();
      for (int iMove = 0; iMove < nMove; ++iMove) {
         attempt = mcMoveManager_[iMove].nAttempt();
         accept  = mcMoveManager_[iMove].nAccept();
         Log::file() << setw(32) << left
              << mcMoveManager_[iMove].className()
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
      UTIL_CHECK(hasWC_);
      UTIL_CHECK(hasMcHamiltonian_);
      UTIL_CHECK(mcState_.isAllocated);
      UTIL_CHECK(!mcState_.hasData);
      
      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (mcState_.w[i].cDField(), system().w().rgrid(i).cDField(), meshSize);
         assignReal<<<nBlocks, nThreads>>>
            (mcState_.wc[i].cDField(), wc_[i].cDField(), meshSize);
      }
      mcState_.mcHamiltonian  = mcHamiltonian_;
      mcState_.mcIdealHamiltonian  = mcIdealHamiltonian_;
      mcState_.mcFieldHamiltonian  = mcFieldHamiltonian_;
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
      
      int nMonomer = system().mixture().nMonomer();
      int meshSize = system().domain().mesh().size();
      
      system().setWRGrid(mcState_.w); 

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      for (int i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (wc_[i].cDField(), mcState_.wc[i].cDField(), meshSize);
      }
      mcHamiltonian_ = mcState_.mcHamiltonian;
      mcIdealHamiltonian_ = mcState_.mcIdealHamiltonian;
      mcFieldHamiltonian_ = mcState_.mcFieldHamiltonian;
      hasMcHamiltonian_ = true;
      hasWC_ = true;
      mcState_.hasData = false;
   }
 
   /*
   * Clear the saved Monte-Carlo state.
   *
   * Used when an attempted Monte-Carlo move is accepted.
   */
   template <int D>
   void McSimulator<D>::clearMcState()
   {  mcState_.hasData = false; }

   /*
   * Compute Monte Carlo Hamiltonian.
   */
   template <int D>
   void McSimulator<D>::computeMcHamiltonian()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().hasCFields());
      UTIL_CHECK(hasWC_);
      hasMcHamiltonian_ = false;

      Mixture<D> const & mixture = system().mixture();
      Domain<D> const & domain = system().domain();

      const int nMonomer = mixture.nMonomer();
      const int meshSize = domain.mesh().size();

      const int np = mixture.nPolymer();
      const int ns = mixture.nSolvent();
      double phi, mu;
      double lnQ = 0.0;

      // Compute polymer ideal gas contributions to lnQ
      if (np > 0) {
         Polymer<D> const * polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-08) {
               lnQ += phi*( -mu + 1.0 )/length;
            }
         }
      }

      // Compute solvent ideal gas contributions to lnQ
      if (ns > 0) {
         Solvent<D> const * solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-8) {
               lnQ += phi*( -mu + 1.0 )/size;
            }
         }
      }
      
      #if 0
      Log::file() << "-lnQ " << -lnQ<< "\n";
      #endif
      
      // lnQ now contains a value per monomer
      
      // Compute field contribution HW
      double HW = 0.0;
      double prefactor;
      // Compute quadratic field contribution to HW
      for (int j = 0; j < nMonomer - 1; ++j) {
         prefactor = -0.5*double(nMonomer)/chiEvals_[j];
         double wSqure = 0;
         wSqure = (double)gpuInnerProduct(wc_[j].cDField(), wc_[j].cDField(), meshSize);
         HW += prefactor * wSqure;
      }
      
      // Subtract average of Langrange multiplier field
      double sum_xi = (double)gpuSum(wc_[nMonomer-1].cDField(), meshSize);
      HW -= sum_xi;
      
      // Normalize HW to equal a value per monomer
      HW /= double(meshSize);
      // Compute final MC Hamiltonian
      mcHamiltonian_ = HW - lnQ;
      const double vSystem  = domain.unitCell().volume();
      const double vMonomer = mixture.vMonomer();
      mcFieldHamiltonian_ = vSystem/vMonomer * HW;
      mcIdealHamiltonian_ = vSystem/vMonomer * lnQ;
      mcHamiltonian_ *= vSystem/vMonomer;
      hasMcHamiltonian_ = true;
      //Log::file()<< "computeMcHamiltonian"<< std::endl;
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

      // Compute chiP = = P*chi*P = P*T
      DMatrix<double> chiP;
      chiP.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            chiP(i, j) = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               chiP(i,j) += P(i,k)*T(k,j);
            }
         }
      }

      // Eigenvalue calculations use data structures and 
      // functions from the Gnu Scientific Library (GSL)

      // Allocate GSL matrix A that will hold a copy of chiP
      gsl_matrix* A = gsl_matrix_alloc(nMonomer, nMonomer);

      // Copy DMatrix<double> chiP to gsl_matrix A
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            gsl_matrix_set(A, i, j, chiP(i, j));
         }
      }

      // Compute eigenvalues and eigenvectors of chiP (or A)
      gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(nMonomer);
      gsl_vector* Avals = gsl_vector_alloc(nMonomer);
      gsl_matrix* Avecs = gsl_matrix_alloc(nMonomer, nMonomer);
      int error;
      error = gsl_eigen_symmv(A, Avals, Avecs, work);

      // Requirements: 
      // - A has exactly one zero eigenvalue, with eigenvector (1,...,1)
      // - All other eigenvalues must be negative.

      // Copy eigenpairs with non-null eigenvalues
      int iNull = -1;  // index for null eigenvalue
      int nNull =  0;  // number of null eigenvalue
      k = 0;           // re-ordered index for non-null eigenvalue
      double val;
      for (i = 0; i < nMonomer; ++i) {
         val = gsl_vector_get(Avals, i);
         if (abs(val) < 1.0E-8) {
            ++nNull;
            iNull = i;
            UTIL_CHECK(nNull <= 1);
         } else {
            chiEvals_[k] = val;
            UTIL_CHECK(val < 0.0);
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
    
      // Set eigenpair with zero eigenvalue 
      i = nMonomer - 1;
      chiEvals_[i] = 0.0;
      for (j = 0; j < nMonomer; ++j) {
         chiEvecs_(i, j) = gsl_matrix_get(Avecs, j, iNull);
      }
      if (chiEvecs_(i, 0) < 0) {
         for (j = 0; j < nMonomer; ++j) {
            chiEvecs_(i, j) = -chiEvecs_(i, j);
         }
      }

      // Normalize all eigenvectors so that the sum of squares = nMonomer
      double vec, norm, prefactor;
      for (i = 0;  i < nMonomer; ++i) {
         norm = 0.0;
         for (j = 0;  j < nMonomer; ++j) {
            vec = chiEvecs_(i, j);
            norm += vec*vec;
         } 
         prefactor = sqrt( double(nMonomer)/norm );
         for (j = 0;  j < nMonomer; ++j) {
            chiEvecs_(i, j) *= prefactor;
         } 
      }

      // Check final eigenvector is (1, ..., 1)
      for (j = 0; j < nMonomer; ++j) {
         UTIL_CHECK(abs(chiEvecs_(nMonomer-1, j) - 1.0) < 1.0E-8);
      }

      #if 0
      // Debugging output
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

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs_ of the projected chi matrix as a basis.
   */
   template <int D>
   void McSimulator<D>::computeWC()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int j, k;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      DArray<RDField<D>> const * currSys = &system().w().rgrid();
      // Loop over eigenvectors (j is an eigenvector index)
      for (j = 0; j < nMonomer; ++j) {
         // Loop over grid points to zero out field wc_[j]
         RDField<D>& wc = wc_[j];
         assignUniformReal<<<nBlocks, nThreads>>>(wc.cDField(), 0, meshSize);
         // Loop over monomer types (k is a monomer index)
         for (k = 0; k < nMonomer; ++k) {
            cudaReal vec;
            vec = (cudaReal)chiEvecs_(j, k)/nMonomer;
            // Loop over grid points
            pointWiseAddScale<<<nBlocks, nThreads>>>
               (wc.cDField(), (*currSys)[k].cDField(), vec, meshSize);
         }
      }
      
      #if 0
      // Debugging output
      std::string filename = "wc";
      system().fieldIo().writeFieldsRGrid(filename, wc_, system().domain().unitCell());
      hasWC_ = true;
      #endif
   }
   
   /*
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void McSimulator<D>::analyzeTrajectory(int min, int max,
                                          std::string classname,
                                          std::string filename)
   {
      // Preconditions
      if (min < 0) UTIL_THROW("min < 0");
      if (max < 0) UTIL_THROW("max < 0");
      if (max < min) UTIL_THROW("max < min!");
      UTIL_CHECK(Analyzer<D>::baseInterval);
      UTIL_CHECK(analyzerManager_.size() > 0);
      
      // Construct TrajectoryReader
      TrajectoryReader<D>* trajectoryReaderPtr;
      trajectoryReaderPtr = trajectoryReaderFactory().factory(classname);
      if (!trajectoryReaderPtr) {
         std::string message;
         message = "Invalid TrajectoryReader class name " + classname;
         UTIL_THROW(message.c_str());
      }
      // Open trajectory file
      Log::file() << "Reading " << filename << std::endl;
      trajectoryReaderPtr->open(filename);
      trajectoryReaderPtr->readHeader();
      // Read Header
      // Main loop over trajectory frames
      Timer timer;
      Log::file() << "Begin main loop" << std::endl;
      bool hasFrame = true;
      timer.start();
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         hasFrame = trajectoryReaderPtr->readFrame();
         if (hasFrame) {
            clearData();
            // Initialize analyzers 
            if (iStep_ == min) analyzerManager_.setup();
            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
               if ((iStep_ % 100) == 0){
                  Log::file() << "Analyzing steps: " << iStep_ << std::endl;
               }
            }
         }
      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;
      int nFrames = iStep_ - min;
      trajectoryReaderPtr->close();
      delete trajectoryReaderPtr;
      // Output results of all analyzers to output files
      analyzerManager_.output();
      // Output time 
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time() 
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames) 
                  << "  sec" << std::endl;
      Log::file() << std::endl;
   }

}
}
#endif
