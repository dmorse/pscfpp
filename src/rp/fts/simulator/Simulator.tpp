#ifndef RP_SIMULATOR_TPP
#define RP_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"

#include <pscf/interaction/Interaction.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/math/IntVec.h>
#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>
#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   Simulator<D,T>::Simulator(typename T::System& system,
                             typename T::Simulator& simulator)
    : hamiltonian_(0.0),
      idealHamiltonian_(0.0),
      fieldHamiltonian_(0.0),
      perturbationHamiltonian_(0.0),
      iStep_(0),
      iTotalStep_(0),
      seed_(0),
      hasHamiltonian_(false),
      hasWc_(false),
      hasCc_(false),
      hasDc_(false),
      systemPtr_(&system),
      simulatorPtr_(&simulator),
      randomPtr_(nullptr),
      vecRandomPtr_(nullptr),
      compressorFactoryPtr_(nullptr),
      compressorPtr_(nullptr),
      perturbationFactoryPtr_(nullptr),
      perturbationPtr_(nullptr),
      rampFactoryPtr_(nullptr),
      rampPtr_(nullptr),
      isAllocated_(false)
   {
      ParamComposite::setClassName("Simulator");
      randomPtr_ = new Random();
      vecRandomPtr_ = new typename T::VecRandom();
      compressorFactoryPtr_ = new typename T::CompressorFactory(system);
      perturbationFactoryPtr_ 
             = new typename T::PerturbationFactory(simulator);
      rampFactoryPtr_ = new typename T::RampFactory(simulator);
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   Simulator<D,T>::~Simulator()
   {
      if (compressorFactoryPtr_) {
         delete compressorFactoryPtr_;
      }
      if (compressorPtr_) {
         delete compressorPtr_;
      }
      if (perturbationFactoryPtr_) {
         delete perturbationFactoryPtr_;
      }
      if (perturbationPtr_) {
         delete perturbationPtr_;
      }
      if (rampFactoryPtr_) {
         delete rampFactoryPtr_;
      }
      if (rampPtr_) {
         delete rampPtr_;
      }
   }

   /*
   * Allocate required memory.
   */
   template <int D, class T>
   void Simulator<D,T>::allocate()
   {
      UTIL_CHECK(!isAllocated_);

      const int nMonomer = system().mixture().nMonomer();

      // Allocate projected chi matrix chiP_ and associated arrays
      chiP_.allocate(nMonomer, nMonomer);
      chiEvals_.allocate(nMonomer);
      chiEvecs_.allocate(nMonomer, nMonomer);
      sc_.allocate(nMonomer);

      // Allocate memory for eignevector components of w and c fields
      wc_.allocate(nMonomer);
      cc_.allocate(nMonomer);
      const IntVec<D> dimensions = system().domain().mesh().dimensions();
      for (int i = 0; i < nMonomer; ++i) {
         wc_[i].allocate(dimensions);
         cc_[i].allocate(dimensions);
      }

      // Allocate memory for components of d (functional derivative)
      dc_.allocate(nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(dimensions);
      }

      // Allocate memory for r-field used as workspace
      tmpField_.allocate(dimensions);

      // Allocate state, if necessary.
      if (!state().isAllocated) {
         state().allocate(nMonomer, dimensions);
      }

      isAllocated_ = true;
   }

   /*
   * Read parameter block.
   *
   * Virtual function - this default version is only used for unit tests.
   */
   template <int D, class T>
   void Simulator<D,T>::readParameters(std::istream &in)
   {
      // Optionally read a random number generator seed
      readRandomSeed(in);

      bool isEnd = false;

      // Optionally read a Compressor block
      readCompressor(in, isEnd);

      // Optionally read a Perturbation block
      readPerturbation(in, isEnd);

      // Optionally read a Ramp block
      readRamp(in, isEnd);
   }

   /*
   * Perform a field theoretic simulation (unimplemented default).
   */
   template <int D, class T>
   void Simulator<D,T>::simulate(int nStep)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D,T>::simulate"); }

   /*
   * Open, read and analyze a trajectory file (unimplemented default).
   */
   template <int D, class T>
   void Simulator<D,T>::analyze(int min, int max,
                                std::string classname,
                                std::string filename)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D,T>::analyze"); }

   /*
   * Clear all local state data (field eigen-components and Hamiltonian)
   */
   template <int D, class T>
   void Simulator<D,T>::clearData()
   {
      hasHamiltonian_ = false;
      hasWc_ = false;
      hasCc_ = false;
      hasDc_ = false;
   }

   /*
   * Compute field theoretic Hamiltonian and its components.
   */
   template <int D, class T>
   void Simulator<D,T>::computeHamiltonian()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().c().hasData());
      UTIL_CHECK(hasWc_);
      hasHamiltonian_ = false;

      Mixture<D,T> const & mixture = system().mixture();
      Domain<D,T> const & domain = system().domain();

      const int nMonomer = mixture.nMonomer();
      const int meshSize = domain.mesh().size();

      const int np = mixture.nPolymer();
      const int ns = mixture.nSolvent();
      double phi, mu;

      // Compute field Hamiltonian contribution
      idealHamiltonian_ = 0.0;

      // Polymer ideal gas contribution
      if (np > 0) {
         Polymer<D,T> const * polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            if (PolymerModel::isThread()) {
               length = polymerPtr->length();
            } else {
               length = (double)polymerPtr->nBead();
            }
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-08) {
               idealHamiltonian_ += phi * ( mu - 1.0 ) / length;
            }
         }
      }

      // Compute solvent ideal gas contributions to idealHamiltonian_
      if (ns > 0) {
         Solvent<D,T> const * solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-8) {
               idealHamiltonian_ += phi * ( mu - 1.0 ) / size;
            }
         }
      }

      // Subtract spatial average of pressure field from idealHamiltonian_
      double sum_xi = Reduce::sum(wc_[nMonomer-1]);
      idealHamiltonian_ -= sum_xi / double(meshSize);

      // idealHamiltonian_  now contains a final value per monomer

      // Compute field Hamiltonian contribution
      fieldHamiltonian_ = 0.0;

      // Quadratic field contribution 
      double prefactor, wSquare;
      for (int j = 0; j < nMonomer - 1; ++j) {
         UTIL_CHECK(tmpField_.capacity() == meshSize);
         // Subtract constant shift sc_[j]
         VecOp::subVS(tmpField_, wc_[j], sc_[j]);
         wSquare = Reduce::sumSq(tmpField_);
         prefactor = -0.5*double(nMonomer)/chiEvals_[j];
         fieldHamiltonian_ += prefactor * wSquare;
      }
      fieldHamiltonian_ /= double(meshSize);

      // Add constant term K/2 per monomer ( K = s = e^{T}chi e/M^2 )
      fieldHamiltonian_ += 0.5*sc_[nMonomer - 1];

      // fieldHamiltonian_ now contains a final value per monomer

      // Compute nMonomerSystem = number of monomers in the system 
      const double vSystem  = domain.unitCell().volume();
      const double vMonomer = mixture.vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;

      // Convert per monomer values to extensive values
      fieldHamiltonian_ *= nMonomerSystem;
      idealHamiltonian_ *= nMonomerSystem;

      hamiltonian_ = idealHamiltonian_ + fieldHamiltonian_;

      // Add perturbationHamiltonian_, if any 
      if (hasPerturbation()) {
         perturbationHamiltonian_ 
            = perturbation().hamiltonian(hamiltonian_);
         hamiltonian_ += perturbationHamiltonian_;
      } else {
         perturbationHamiltonian_ = 0.0;
      }

      hasHamiltonian_ = true;
   }

   /*
   * Compute eigenvalues and eigenvectors of the projected chi matrix.
   */
   template <int D, class T>
   void Simulator<D,T>::analyzeChi()
   {
      UTIL_CHECK(isAllocated_);

      const int nMonomer = system().mixture().nMonomer();
      DMatrix<double> const & chi = system().interaction().chi();
      double d = 1.0/double(nMonomer);
      int i, j, k;

      // Compute orthogonal projection matrix P
      DMatrix<double> P;
      P.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            P(i,j) = -d;
         }
         P(i,i) += 1.0;
      }

      // Compute CP = chi*P (temporary matrix)
      DMatrix<double> CP;
      CP.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            CP(i, j) = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               CP(i,j) += chi(i,k)*P(k,j);
            }
         }
      }

      // Compute chiP = = P*chi*P = P*CP
      DMatrix<double> chiP;
      chiP.allocate(nMonomer, nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            chiP(i, j) = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               chiP(i,j) += P(i,k)*CP(k,j);
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
      UTIL_CHECK(error == 0);

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
         if (std::abs(val) < 1.0E-8) {
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
         UTIL_CHECK(std::abs(chiEvecs_(nMonomer-1, j) - 1.0) < 1.0E-8);
      }

      // Compute vector s in monomer basis
      DArray<double> s;
      s.allocate(nMonomer);
      for (i = 0; i < nMonomer; ++i) {
         s[i] = 0.0;
         for (j = 0; j < nMonomer; ++j) {
           s[i] += chi(i,j);
         }
         s[i] = s[i]/double(nMonomer);
      }

      // Compute components of s in eigenvector basis -> sc_
      for (i = 0; i < nMonomer; ++i) {
         sc_[i] = 0.0;
         for (j = 0; j < nMonomer; ++j) {
            sc_[i] += chiEvecs_(i,j)*s[j];
         }
         sc_[i] = sc_[i]/double(nMonomer);
      }

   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs_ of the projected chi matrix as a basis.
   */
   template <int D, class T>
   void Simulator<D,T>::computeWc()
   {
      UTIL_CHECK(isAllocated_);

      double vec;
      int i, j;
      const int nMonomer = system().mixture().nMonomer();

      // Loop over eigenvectors (i is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {

         // Initialize field wc_[i] to zero
         typename T::RField& Wc = wc_[i];
         VecOp::eqS(Wc, 0.0);

         // Loop over monomer types (j is a monomer index)
         for (j = 0; j < nMonomer; ++j) {
            vec = chiEvecs_(i, j)/double(nMonomer);
            VecOp::addEqVc(Wc, system().w().rgrid(j), vec);
         }
      }

      hasWc_ = true;
   }

   /*
   * Compute the eigenvector components of the c-fields, using the
   * eigenvectors chiEvecs_ of the projected chi matrix as a basis.
   */
   template <int D, class T>
   void Simulator<D,T>::computeCc()
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().c().hasData());

      double vec;
      int i, j;
      const int nMonomer = system().mixture().nMonomer();

      // Loop over eigenvectors (i is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {

         // Initialize field cc_[i] to zero
         typename T::RField& Cc = cc_[i];
         VecOp::eqS(Cc, 0.0);

         // Loop over monomer types (j is a monomer index)
         for (j = 0; j < nMonomer; ++j) {
            vec = chiEvecs_(i, j);
            VecOp::addEqVc(Cc, system().c().rgrid(j), vec);
         }
      }

      hasCc_ = true;
   }

   /*
   * Compute d fields, i.e., functional derivatives of H[W].
   */
   template <int D, class T>
   void Simulator<D,T>::computeDc()
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      if (!hasWc_) computeWc();
      if (!hasCc_) computeCc();

      // Local constants and variables
      const int nMonomer = system().mixture().nMonomer();
      const double vMonomer = system().mixture().vMonomer();
      const double a = 1.0/vMonomer;
      double b, s;

      // Loop over composition eigenvectors (exclude the last)
      for (int i = 0; i < nMonomer - 1; ++i) {
         typename T::RField& Dc = dc_[i];
         typename T::RField const & Wc = wc_[i];
         typename T::RField const & Cc = cc_[i];
         b = -1.0*a*double(nMonomer)/chiEvals_[i];
         s = -1.0*b*sc_[i];
         VecOp::addVcVcS(Dc, Cc, a, Wc, b, s);
      }

      // Add derivatives arising from a perturbation (if any).
      if (hasPerturbation()) {
         perturbation().incrementDc(dc_);
      }

      hasDc_ = true;
   }

   /*
   * Save the current state prior to a next move.
   *
   * Invoked before each attempted move.
   */
   template <int D, class T>
   void Simulator<D,T>::saveState()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(hasWc());
      UTIL_CHECK(state().isAllocated);
      UTIL_CHECK(!state().hasData);

      const int nMonomer = system().mixture().nMonomer();

      // Save all w-field components
      for (int i = 0; i < nMonomer; ++i) {
         state().w[i] = system().w().rgrid(i);
         state().wc[i] = wc_[i];
      }

      // Save cc, if needed
      if (state().needsCc) {
         UTIL_CHECK(hasCc());
         UTIL_CHECK(state().cc.isAllocated());
         for (int i = 0; i < nMonomer; ++i) {
            state().cc[i] = cc_[i];
         }
      }

      // Save dc, if needed
      if (state().needsDc) {
         UTIL_CHECK(hasDc());
         UTIL_CHECK(state().dc.isAllocated());
         for (int i = 0; i < nMonomer - 1; ++i) {
            state().dc[i] = dc_[i];
         }
      }

      // Save Hamiltonian and its components, if needed
      if (state().needsHamiltonian){
         UTIL_CHECK(hasHamiltonian());
         state().hamiltonian  = hamiltonian();
         state().idealHamiltonian  = idealHamiltonian();
         state().fieldHamiltonian  = fieldHamiltonian();
         state().perturbationHamiltonian  = perturbationHamiltonian();
      }

      if (hasPerturbation()) {
         perturbation().saveState();
      }

      state().hasData = true;
   }

   /*
   * Restore a saved simulation state.
   *
   * Invoked after the compressor fails to converge or an attempted
   * Monte-Carlo move is rejected.
   */
   template <int D, class T>
   void Simulator<D,T>::restoreState()
   {
      UTIL_CHECK(state().isAllocated);
      UTIL_CHECK(state().hasData);
      const int nMonomer = system().mixture().nMonomer();

      // Restore monomer w-field components
      system().w().setRGrid(state().w);

      // Restore wc_ field components
      for (int i = 0; i < nMonomer; ++i) {
         wc_[i] = state().wc[i];
      }
      hasWc_ = true;

      // Restore cc_ c-field components, if needed
      if (state().needsCc) {
         for (int i = 0; i < nMonomer; ++i) {
            cc_[i] = state().cc[i];
         }
         hasCc_ = true;
      }

      // Restore dc_ fieldc components, if needed
      if (state().needsDc) {
         for (int i = 0; i < nMonomer - 1; ++i) {
            dc_[i] = state().dc[i];
         }
         hasDc_ = true;
      }

      // Restore Hamiltonian and its components, if needed
      if (state().needsHamiltonian){
         hamiltonian_ = state().hamiltonian;
         idealHamiltonian_ = state().idealHamiltonian;
         fieldHamiltonian_ = state().fieldHamiltonian;
         perturbationHamiltonian_ = state().perturbationHamiltonian;
         hasHamiltonian_ = true;
      }

      // Restore internal state of perturbation, if any
      if (hasPerturbation()) {
         perturbation().restoreState();
      }

      // Clear internal state of SimState object
      state().hasData = false;
   }

   /*
   * Clear the saved system state.
   *
   * Invoked when an attempted move is accepted.
   */
   template <int D, class T>
   void Simulator<D,T>::clearState()
   {  state().hasData = false; }

   /*
   * Output all timer results.
   */
   template <int D, class T>
   void Simulator<D,T>::outputTimers(std::ostream& out) const
   {
      UTIL_CHECK(compressorPtr_);
      outputMdeCounter(out);
      compressorPtr_->outputTimers(out);
   }

   /*
   * Output modified diffusion equation (MDE) counter.
   */
   template <int D, class T>
   void Simulator<D,T>::outputMdeCounter(std::ostream& out) const
   {
      UTIL_CHECK(compressorPtr_);
      out << "MDE counter   "
          << compressorPtr_->mdeCounter() << std::endl;
      out << std::endl;
   }

   /*
   * Clear all timers.
   */
   template <int D, class T>
   void Simulator<D,T>::clearTimers()
   {
      UTIL_CHECK(hasCompressor());
      compressor().clearTimers();
   }

   // Protected Functions

   /*
   * Optionally read RNG seed, initialize random number generators.
   */
   template <int D, class T>
   void Simulator<D,T>::readRandomSeed(std::istream& in)
   {
      // Optionally read a random number generator seed
      seed_ = 0;
      readOptional(in, "seed", seed_);

      // Set random number generator seeds on CPU and GPU
      // Default value seed_ = 0 uses the clock time.
      random().setSeed(seed_);
      initializeVecRandom();
   }

   // Functions related to a Compressor

   /*
   * Optionally read a Compressor parameter file block.
   */
   template <int D, class T>
   void Simulator<D,T>::readCompressor(std::istream& in, bool& isEnd)
   {
      if (!isEnd) {
         UTIL_CHECK(compressorFactoryPtr_);
         UTIL_CHECK(!hasCompressor());
         std::string className;
         compressorPtr_ =
            compressorFactory().readObjectOptional(in, *this,
                                                   className, isEnd);
      }
      if (!compressorPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Compressor{ [absent] }\n";
      }
   }

   /*
   * Get the Compressor factory by reference.
   */
   template <int D, class T>
   typename T::CompressorFactory& Simulator<D,T>::compressorFactory()
   {
      UTIL_CHECK(compressorFactoryPtr_);
      return *compressorFactoryPtr_;
   }

   // Functions related to a Perturbation

   /*
   * Optionally read a Perturbation parameter file block.
   */
   template <int D, class T>
   void Simulator<D,T>::readPerturbation(std::istream& in, bool& isEnd)
   {
      if (!isEnd) {
         UTIL_CHECK(perturbationFactoryPtr_);
         UTIL_CHECK(!hasPerturbation());
         std::string className;
         perturbationPtr_ =
            perturbationFactory().readObjectOptional(in, *this,
                                                     className, isEnd);
      }
      if (!perturbationPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Perturbation{ [absent] }\n";
      }
   }

   /*
   * Get the Perturbation factory by reference.
   */
   template <int D, class T>
   typename T::PerturbationFactory& Simulator<D,T>::perturbationFactory()
   {
      UTIL_CHECK(perturbationFactoryPtr_);
      return *perturbationFactoryPtr_;
   }

   /*
   * Set the associated Perturbation object.
   */
   template <int D, class T>
   void Simulator<D,T>::setPerturbation(typename T::Perturbation* ptr)
   {
      UTIL_CHECK(ptr);
      perturbationPtr_ = ptr;
   }

   // Functions related to a Ramp

   /*
   * Optionally read a Ramp parameter file block.
   */
   template <int D, class T>
   void Simulator<D,T>::readRamp(std::istream& in, bool& isEnd)
   {
      if (!isEnd) {
         UTIL_CHECK(rampFactoryPtr_);
         UTIL_CHECK(!hasRamp());
         std::string className;
         rampPtr_ =
            rampFactory().readObjectOptional(in, *this, className, isEnd);
      }
      if (!rampPtr_ && ParamComponent::echo()) {
         Log::file() << indent() << "  Ramp{ [absent] }\n";
      }
   }

   /*
   * Get the Ramp factory by reference.
   */
   template <int D, class T>
   typename T::RampFactory& Simulator<D,T>::rampFactory()
   {
      UTIL_CHECK(rampFactoryPtr_);
      return *rampFactoryPtr_;
   }

   /*
   * Set the associated Ramp object.
   */
   template <int D, class T>
   void Simulator<D,T>::setRamp(typename T::Ramp* ptr)
   {
      UTIL_CHECK(ptr);
      rampPtr_ = ptr;
   }

} // namespace Rp
} // namespace Pscf
#endif
