#ifndef CPC_SIMULATOR_TPP
#define CPC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Contains declarations needed in some other templates
#include <pscf/cpu/complex.h>

#include "Simulator.h"
#include <cpc/system/System.h>
#include <cpc/solvers/Mixture.h>
#include <cpc/solvers/Polymer.h>
#include <cpc/solvers/Solvent.h>
#include <cpc/field/Domain.h>
#include <cpc/fts/step/Step.h>
#include <cpc/fts/step/StepFactory.h>
#if 0
#include <cpc/fts/analyzer/AnalyzerFactory.h>
#include <cpc/fts/trajectory/TrajectoryReader.h>
#include <cpc/fts/trajectory/TrajectoryReaderFactory.h>
#endif
#include <pscf/interaction/Interaction.h>
#include <util/random/Random.h>
#include <util/misc/Timer.h>
#include <util/global.h>

// Gnu scientific library
#include <gsl/gsl_eigen.h>

#include <complex>

namespace Pscf {
namespace Cpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Simulator<D>::Simulator(System<D>& system)
    : //analyzerManager_(*this, system),
      random_(),
      wc_(),
      cc_(),
      dc_(),
      u_(),
      evecs_(),
      evals_(),
      isPositiveEval_(),
      hamiltonian_(0.0),
      idealHamiltonian_(0.0),
      fieldHamiltonian_(0.0),
      systemPtr_(&system),
      stepPtr_(nullptr),
      stepFactoryPtr_(nullptr),
      //trajectoryReaderFactoryPtr_(nullptr)
      iStep_(0),
      iTotalStep_(0), 
      seed_(0),
      hasHamiltonian_(false),
      hasWc_(false),
      hasCc_(false),
      hasDc_(false),
      isAllocated_(false)
   {  
      ParamComposite::setClassName("Simulator"); 
      stepFactoryPtr_ = new StepFactory<D>(*this);
      // trajectoryReaderFactoryPtr_
      //       = new TrajectoryReaderFactory<D>(system);
   }

   /*
   * Destructor.
   */
   template <int D>
   Simulator<D>::~Simulator()
   {
      if (stepPtr_) {
         delete stepPtr_;
      }
      if (stepFactoryPtr_) {
         delete stepFactoryPtr_;
      }
      #if 0
      if (trajectoryReaderFactoryPtr_) {
         delete trajectoryReaderFactoryPtr_;
      }
      #endif
   }

   /*
   * Allocate required memory.
   */
   template <int D>
   void Simulator<D>::allocate()
   {
      UTIL_CHECK(!isAllocated_);

      const int nMonomer = system().mixture().nMonomer();

      // Allocate interaction matrix u_ and associated arrays
      u_.allocate(nMonomer, nMonomer);
      evecs_.allocate(nMonomer, nMonomer);
      evals_.allocate(nMonomer);
      isPositiveEval_.allocate(nMonomer);

      // Allocate memory for eignevector components of w and c fields
      wc_.allocate(nMonomer);
      cc_.allocate(nMonomer);
      dc_.allocate(nMonomer);
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      for (int i = 0; i < nMonomer; ++i) {
         wc_[i].allocate(dimensions);
         cc_[i].allocate(dimensions);
         dc_[i].allocate(dimensions);
      }

      isAllocated_ = true;
   }
   
   /*
   * Read parameter file block for a BD simulator.
   */
   template <int D>
   void Simulator<D>::readParameters(std::istream &in)
   {
      // Optionally read a random seed value
      readRandomSeed(in);
  
      // Optionally read a Step block
      bool isEnd = false;
      std::string className;
      stepPtr_ =
         stepFactoryPtr_->readObjectOptional(in, *this, className, 
                                             isEnd);
      if (!hasStep() && ParamComponent::echo()) {
         Log::file() << indent() << "  Step{ [absent] }\n";
      }

      #if 0 
      // Optionally read an AnalyzerManager
      Analyzer<D>::baseInterval = 0; // default value
      readParamCompositeOptional(in, analyzerManager_);
      #endif

   }

   /*
   * Setup before main loop of a BD simulation.
   */
   template <int D>
   void Simulator<D>::setup(int nStep)
   {
      UTIL_CHECK(system().w().hasData());

      // Eigenanalysis of the interaction matrix
      analyzeInteraction();

      // Solve MDE and compute c-fields for the intial state
      system().compute();

      // Compute field components and Hamiltonian for initial state.
      computeWc();
      computeCc();
      computeDc();
      computeHamiltonian();

      step().setup();
      #if 0
      if (analyzerManager_.size() > 0){
         analyzerManager_.setup();
      }
      #endif

   }

   /*
   * Perform a field theoretic simulation (unimplemented by base class).
   */
   template <int D>
   void Simulator<D>::simulate(int nStep)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D>::simulate"); }

   #if 0
   /*
   * Perform a field theoretic MC simulation of nStep steps.
   */
   template <int D>
   void Simulator<D>::simulate(int nStep)
   {
      UTIL_CHECK(hasStep());
      UTIL_CHECK(system().w().hasData());

      // Initial setup
      setup(nStep);

      // Main simulation loop
      Timer timer;
      Timer analyzerTimer;
      timer.start();
      iStep_ = 0;

      #if 0
      // Analysis for initial state (if any)
      analyzerTimer.start();
      if (analyzerManager_.size() > 0){
         analyzerManager_.sample(iStep_);
      }
      analyzerTimer.stop();
      #endif

      for (iTotalStep_ = 0; iTotalStep_ < nStep; ++iTotalStep_) {

         // Take a step (modifies W fields)
         step().step();
         iStep_++;

         #if 0
         // Analysis (if any)
         analyzerTimer.start();
         if (Analyzer<D>::baseInterval != 0) {
            if (analyzerManager_.size() > 0) {
               if (iStep_ % Analyzer<D>::baseInterval == 0) {
                  analyzerManager_.sample(iStep_);
               }
            }
         }
         analyzerTimer.stop();
         #endif

      }

      timer.stop();
      double time = timer.time();
      double analyzerTime = analyzerTimer.time();

      #if 0
      // Output results analyzers to files
      if (Analyzer<D>::baseInterval > 0){
         analyzerManager_.output();
      }
      #endif

      // Output times for the simulation run
      Log::file() << std::endl;
      Log::file() << "nStep               " << nStep << std::endl;
      Log::file() << "Total run time      " << time
                  << " sec" << std::endl;
      double rStep = double(nStep);
      Log::file() << "time / nStep        " <<  time / rStep
                  << " sec" << std::endl;
      Log::file() << "Analyzer run time   " << analyzerTime
                  << " sec" << std::endl;
      Log::file() << std::endl;

   }
   #endif

   #if 0
   /*
   * Open, read and analyze a trajectory file (unimplemented by base class).
   */
   template <int D>
   void Simulator<D>::analyze(int min, int max,
                              std::string classname,
                              std::string filename)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D>::analyze"); }
   #endif

   #if 0
   /*
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void Simulator<D>::analyze(int min, int max,
                                std::string classname,
                                std::string filename)
   {
      // Preconditions
      UTIL_CHECK(min >= 0);
      UTIL_CHECK(max >= min);
      UTIL_CHECK(Analyzer<D>::baseInterval > 0);
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
      trajectoryReaderPtr->open(filename);
      trajectoryReaderPtr->readHeader();

      // Main loop over trajectory frames
      Timer timer;
      bool hasFrame;
      timer.start();
      hasFrame = trajectoryReaderPtr->readFrame();
      
      for (iStep_ = 0; iStep_ <= max && hasFrame; ++iStep_) {
         if (hasFrame) {
            clearData();

            // Initialize analyzers
            if (iStep_ == min) {
               //analyzerManager_.setup();
               setup(iStep_);
            }

            // Sample property values only for iStep >= min
            if (iStep_ >= min) {
               analyzerManager_.sample(iStep_);
            }
         }
         
         hasFrame = trajectoryReaderPtr->readFrame();
      }
      timer.stop();
      Log::file() << "end main loop" << std::endl;
      int nFrames = iStep_ - min;
      trajectoryReaderPtr->close();
      delete trajectoryReaderPtr;

      // Output results of all analyzers to output files
      analyzerManager_.output();

      // Output number of frames and times
      Log::file() << std::endl;
      Log::file() << "# of frames   " << nFrames << std::endl;
      Log::file() << "run time      " << timer.time()
                  << "  sec" << std::endl;
      Log::file() << "time / frame " << timer.time()/double(nFrames)
                  << "  sec" << std::endl;
      Log::file() << std::endl;

   }
   #endif

   /*
   * Clear all state data (Hamiltonian, eigen-components of w, c, and d)
   */
   template <int D>
   void Simulator<D>::clearData()
   {
      hasHamiltonian_ = false;
      hasWc_ = false;
      hasCc_ = false;
      hasDc_ = false;
   }

   /*
   * Compute field theoretic Hamiltonian H[W].
   */
   template <int D>
   void Simulator<D>::computeHamiltonian()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().c().hasData());
      UTIL_CHECK(hasWc_);
      hasHamiltonian_ = false;

      Mixture<D> const & mixture = system().mixture();
      Domain<D> const & domain = system().domain();

      const int nMonomer = mixture.nMonomer();
      const int meshSize = domain.mesh().size();

      const int np = mixture.nPolymer();
      const int ns = mixture.nSolvent();

      // Compute number of monomers in the system (nMonomerSystem)
      const double vSystem  = domain.unitCell().volume();
      const double vMonomer = mixture.vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;

      // Compute polymer ideal gas contributions to lnQ
      std::complex<double> phi, mu, lnQ;
      lnQ = std::complex<double>(0.0, 0.0);
      if (np > 0) {
         Polymer<D> const * polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture.polymer(i);
            phi = polymerPtr->phi();
            if (std::abs(phi) > 1.0E-08) {
               mu = polymerPtr->mu();
               if (PolymerModel::isThread()) {
                  length = polymerPtr->length();
               } else {
                  length = (double)polymerPtr->nBead();
               }
               lnQ += phi*( -mu + 1.0 )/length;
               // Recall: mu = ln(phi/q)
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
            if (std::abs(phi) > 1.0E-8) {
               mu = solventPtr->mu();
               size = solventPtr->size();
               lnQ += phi*( -mu + 1.0 )/size;
               // Recall: mu = ln(phi/q)
            }
         }
      }
      idealHamiltonian_ = -1.0 * nMonomerSystem * lnQ;

      // Compute quadratic field contribution to fieldHamiltonian_
      fieldHamiltonian_ = std::complex<double>(0.0, 0.0);
      std::complex<double> w;
      double prefactor;
      int i, j;
      for (j = 0; j < nMonomer; ++j) {
         CField<D,CPT> const & Wc = wc_[j];
         prefactor = -0.5*double(nMonomer)/evals_[j];
         for (i = 0; i < meshSize; ++i) {
            assign(w, Wc[i]);
            fieldHamiltonian_ += prefactor * w * w;
         }
      }
      fieldHamiltonian_ /= double(meshSize);
      fieldHamiltonian_ *= nMonomerSystem;

      hamiltonian_ = idealHamiltonian_ + fieldHamiltonian_;

      hasHamiltonian_ = true;
   }

   /*
   * Construct and diagonalize the interaction matrix, U.
   */ 
   template <int D>
   void Simulator<D>::analyzeInteraction()
   {
      UTIL_CHECK(isAllocated_);

      const int nMonomer = system().mixture().nMonomer();
      DMatrix<double> const & chi = system().interaction().chi();
      double const & zeta = system().interaction().zeta();
      int i, j;

      // Compute u_ 
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            u_(i, j) = chi(i,j) + zeta;
         }
      }

      // Eigenvalue calculations use data structures and
      // functions from the Gnu Scientific Library (GSL)

      // Allocate GSL matrix A that will hold a copy of U matrix
      gsl_matrix* A = gsl_matrix_alloc(nMonomer, nMonomer);

      // Copy elements of DMatrix<double> u_ to gsl_matrix A
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            gsl_matrix_set(A, i, j, u_(i, j));
         }
      }

      // Compute eigenvalues and eigenvectors of u_ (or A)
      gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(nMonomer);
      gsl_vector* Avals = gsl_vector_alloc(nMonomer);
      gsl_matrix* Avecs = gsl_matrix_alloc(nMonomer, nMonomer);
      int error;
      error = gsl_eigen_symmv(A, Avals, Avecs, work);
      UTIL_CHECK(error == 0);

      // Copy eigenvalues and eigenvectors
      double val;
      for (i = 0; i < nMonomer; ++i) {
         val = gsl_vector_get(Avals, i);
         evals_[i] = val;
         if (val > 0.0) {
            isPositiveEval_[i] = true;
         } else {
            isPositiveEval_[i] = false;
         }
         for (j = 0; j < nMonomer; ++j) {
            evecs_(i, j) = gsl_matrix_get(Avecs, j, i);
         }
         if (evecs_(i, 0) < 0.0) {
            for (j = 0; j < nMonomer; ++j) {
               evecs_(i, j) = -evecs_(i, j);
            }
         }
      }

      // Normalize all eigenvectors so that the sum of squares = nMonomer
      double vec, norm, prefactor;
      for (i = 0;  i < nMonomer; ++i) {
         norm = 0.0;
         for (j = 0;  j < nMonomer; ++j) {
            vec = evecs_(i, j);
            norm += vec*vec;
         }
         prefactor = sqrt( double(nMonomer)/norm );
         for (j = 0;  j < nMonomer; ++j) {
            evecs_(i, j) *= prefactor;
         }
      }

      #if 0
      // Debugging output
      for (i = 0; i < nMonomer; ++i) {
         Log::file() << "Eigenpair " << i << "\n";
         Log::file() << "value  =  " << evals_[i] << "\n";
         Log::file() << "vector = [ ";
         for (j = 0; j < nMonomer; ++j) {
            Log::file() << evecs_(i, j) << "   ";
         }
         Log::file() << "]\n";
         Log::file() << " sc[i] = " << sc_[i] << std::endl;
      }
      #endif

   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors evecs_ of the interaction matrix as a basis.
   */
   template <int D>
   void Simulator<D>::computeWc()
   {
      UTIL_CHECK(isAllocated_);

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      fftw_complex product;
      double vec;
      int i, j, k;

      // Loop over eigenvectors (j is an eigenvector index)
      for (j = 0; j < nMonomer; ++j) {

         // Loop over grid points to zero out field wc_[j]
         CField<D,CPT>& Wc = wc_[j];
         for (i = 0; i < meshSize; ++i) {
            assign(Wc[i], 0.0);
         }

         // Loop over monomer types (k is a monomer index)
         for (k = 0; k < nMonomer; ++k) {
            vec = evecs_(j, k)/double(nMonomer);

            // Loop over grid points
            CField<D,CPT> const & Wr = system().w().field(k);
            for (i = 0; i < meshSize; ++i) {
               mul(product, Wr[i], vec);
               addEq(Wc[i], product);
               // Wc[i] += vec*Wr[i];
            }

         }
      }

      hasWc_ = true;
   }

   /*
   * Compute the eigenvector components of the c-fields, using the
   * eigenvectors evecs_ of the interaction matrix as a basis.
   */
   template <int D>
   void Simulator<D>::computeCc()
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().c().hasData());

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();

      // Loop over eigenvectors (i is an eigenvector index)
      fftw_complex product;
      double vec;
      int i, j, k;
      for (i = 0; i < nMonomer; ++i) {

         // Initialize field cc_[i] to zero
         CField<D,CPT>& Cc = cc_[i];
         for (k = 0; k < meshSize; ++k) {
            assign(Cc[k], 0.0);
         }

         // Loop over monomer types
         for (j = 0; j < nMonomer; ++j) {
            CField<D,CPT> const & Cr = system().c().field(j);
            vec = evecs_(i, j);

            // Loop over grid points
            for (k = 0; k < meshSize; ++k) {
               //Cc[k] += Cr[k]*vec;
               mul(product, Cr[k], vec);
               addEq(Cc[k], product);
            }

         }
      }

      hasCc_ = true;
   }

   /*
   * Compute d fields, i.e., functional derivatives of H[W].
   */
   template <int D>
   void Simulator<D>::computeDc()
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      if (!hasWc_) computeWc();
      if (!hasCc_) computeCc();

      // Local constants and variables
      const int meshSize = system().domain().mesh().size();
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      const double a = 1.0/vMonomer;
      fftw_complex z;
      double b;
      int i, k;

      // Compute derivatives for standard Hamiltonian
      // Loop over composition eigenvectors (exclude the last)
      for (i = 0; i < nMonomer - 1; ++i) {
         CField<D,CPT>& Dc = dc_[i];
         CField<D,CPT> const & Wc = wc_[i];
         CField<D,CPT> const & Cc = cc_[i]; 
         b = -1.0*double(nMonomer)/evals_[i];
         // Loop over grid points
         for (k = 0; k < meshSize; ++k) {
            //Dc[k] = a*( b*Wc[k] + Cc[k] );
            mul(z, Wc[k], b);
            addEq(z, Cc[k]);
            mul(Dc[k], z, a);
         }
      }

      hasDc_ = true;
   }
   
   // Private Functions

   /*
   * Optionally read a random number generator seed.
   */
   template<int D>
   void Simulator<D>::readRandomSeed(std::istream& in)
   {
      // Optionally read a random number generator seed
      seed_ = 0;
      readOptional(in, "seed", seed_);

      // Set random number generator seed
      // Default value seed_ = 0 uses the clock time.
      random().setSeed(seed_);
   }

} // Cpc
} // Pscf
#endif
