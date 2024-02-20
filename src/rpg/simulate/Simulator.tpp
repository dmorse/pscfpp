#ifndef PSPG_SIMULATOR_TPP
#define PSPG_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"
#include <rpg/System.h>
#include <rpg/compressor/Compressor.h>
#include <util/misc/Timer.h>
#include <util/random/Random.h>
#include <util/global.h>
#include <gsl/gsl_eigen.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Simulator<D>::Simulator(System<D>& system)
   : random_(),
     cudaRandom_(),
     iStep_(0),
     hasHamiltonian_(false),
     hasWc_(false),
     hasCc_(false),
     systemPtr_(&system),
     isAllocated_(false)
   { setClassName("Simulator"); }

   /*
   * Destructor.
   */
   template <int D>
   Simulator<D>::~Simulator()
   {}

   /* 
   * Allocate required memory.
   */
   template <int D>
   void Simulator<D>::allocate()
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
      const int meshSize = system().domain().mesh().size();
      for (int i = 0; i < nMonomer; ++i) {
         wc_[i].allocate(meshSize);
         cc_[i].allocate(meshSize);
      }

      // Allocate memory for components of d (functional derivative)
      dc_.allocate(nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(meshSize);
      }
      
      // Allocate memory for single eignevector components of w after constant shift
      wcs_.allocate(meshSize);
      
      isAllocated_ = true;
   }

   /* 
   * Virtual function to read parameters - unimplemented.
   */
   template <int D>
   void Simulator<D>::readParameters(std::istream &in)
   {  UTIL_THROW("Error: Unimplemented Simulator<D>::readParameters"); } 

   /*
   * Perform a field theoretic simulation of nStep steps.
   */
   template <int D>
   void Simulator<D>::simulate(int nStep)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D>::simulate"); }

   /*
   * Open, read and analyze a trajectory file
   */
   template <int D>
   void Simulator<D>::analyze(int min, int max,
                              std::string classname,
                              std::string filename)
   {  UTIL_THROW("Error: Unimplemented function Simulator<D>::analyze"); }
   
   /*
   * Compute field theoretic Hamiltonian H[W].
   */
   template <int D>
   void Simulator<D>::computeHamiltonian()
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().hasCFields());
      UTIL_CHECK(hasWc_);
      hasHamiltonian_ = false;

      Mixture<D> const & mixture = system().mixture();
      Domain<D> const & domain = system().domain();

      const int nMonomer = mixture.nMonomer();
      const int meshSize = domain.mesh().size();
      
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      const int np = mixture.nPolymer();
      const int ns = mixture.nSolvent();
      double phi, mu;

      // Compute polymer ideal gas contributions to lnQ
      double lnQ = 0.0;
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
      
      // Add average of pressure field wc_[nMonomer-1] to lnQ
      double sum_xi = 
           (double)gpuSum(wc_[nMonomer-1].cField(), meshSize);
      lnQ += sum_xi/double(meshSize);
      
      // lnQ now contains a value per monomer

      // Initialize field contribution HW

      // Compute quadratic field contribution to HW
      double HW = 0.0;
      double prefactor, s;
      int j;
      for (j = 0; j < nMonomer - 1; ++j) {
         prefactor = -0.5*double(nMonomer)/chiEvals_[j];
         // Obtain constant shift
         s = sc_[j];
         // Subtract of constat shift s
         assignRealSubtractDouble<<<nBlocks, nThreads>>>
            (wcs_.cField(), wc_[j].cField(), s, meshSize);
         // Compute quadratic field contribution to HW
         double wSqure = 0;
         wSqure = 
              (double)gpuInnerProduct(wcs_.cField(), 
                                      wcs_.cField(), meshSize);
         HW += prefactor * wSqure;
      }
      
      // Normalize HW to equal a value per monomer
      HW /= double(meshSize);

      // Add constant term K/2 per monomer (K=s=e^{T}chi e/M^2)
      HW += 0.5*sc_[nMonomer - 1];

      // Compute number of monomers in the system (nMonomerSystem)      
      const double vSystem  = domain.unitCell().volume();
      const double vMonomer = mixture.vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;

      // Compute final Hamiltonian components
      fieldHamiltonian_ = nMonomerSystem * HW;
      idealHamiltonian_ = -1.0 * nMonomerSystem * lnQ;
      hamiltonian_ = idealHamiltonian_ + fieldHamiltonian_;

      hasHamiltonian_ = true;
   }

   template <int D>
   void Simulator<D>::analyzeChi()
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
         Log::file() << " sc[i] = " << sc_{i] << std::endl;
      }
      #endif

   }

   /*
   * Compute the eigenvector components of the w fields, using the
   * eigenvectors chiEvecs_ of the projected chi matrix as a basis.
   */
   template <int D>
   void Simulator<D>::computeWc()
   {
      UTIL_CHECK(isAllocated_);

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i,j;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      DArray<RField<D>> const * Wr = &system().w().rgrid();
      // Loop over eigenvectors (j is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {

         // Loop over grid points to zero out field wc_[j]
         RField<D>& Wc = wc_[i];
         assignUniformReal<<<nBlocks, nThreads>>>(Wc.cField(), 0, meshSize);

         // Loop over monomer types (k is a monomer index)
         for (j = 0; j < nMonomer; ++j) {
            cudaReal vec;
            vec = (cudaReal) chiEvecs_(i, j)/nMonomer;

            // Loop over grid points
            pointWiseAddScale<<<nBlocks, nThreads>>>
               (Wc.cField(), (*Wr)[j].cField(), vec, meshSize);
         }
      }
      
      #if 0
      // Debugging output
      std::string filename = "wc";
      system().fieldIo().writeFieldsRGrid(filename, wc_, 
                                          system().domain().unitCell());
      #endif

      hasWc_ = true;
   }
   
   /*
   * Compute the eigenvector components of the c-fields, using the
   * eigenvectors chiEvecs_ of the projected chi matrix as a basis.
   */
   template <int D>
   void Simulator<D>::computeCc()
   {
      // Preconditions
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().hasCFields());

      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);

      DArray<RField<D>> const * Cr = &system().c().rgrid();
      // Loop over eigenvectors (i is an eigenvector index)
      for (i = 0; i < nMonomer; ++i) {
         
         // Set cc_[i] to zero
         RField<D>& Cc = cc_[i];
         assignUniformReal<<<nBlocks, nThreads>>>(Cc.cField(), 0, meshSize);

         // Loop over monomer types 
         for (j = 0; j < nMonomer; ++j) {
            cudaReal vec;
            vec = (cudaReal)chiEvecs_(i, j);
            
            // Loop over grid points
            pointWiseAddScale<<<nBlocks, nThreads>>>
               (Cc.cField(), (*Cr)[j].cField(), vec, meshSize);
         }
      }
      
      #if 0
      // Debugging output
      std::string filename = "cc";
      system().fieldIo().writeFieldsRGrid(filename, cc_, 
                                          system().domain().unitCell());
      #endif

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
      double b, s;
      int i;
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Loop over composition eigenvectors (exclude the last)
      for (i = 0; i < nMonomer - 1; ++i) {
         RField<D>& Dc = dc_[i];
         RField<D> const & Wc = wc_[i];
         RField<D> const & Cc = cc_[i];
         b = -1.0*double(nMonomer)/chiEvals_[i];
         s = sc_[i];
          
         // Loop over grid points
         computeDField<<<nBlocks, nThreads>>>
            (Dc.cField(), Wc.cField(), Cc.cField(), a, b, s, meshSize);
      }

      hasDc_ = true;
   }

   /*
   * Output all timer results.
   */ 
   template<int D>
   void Simulator<D>::outputTimers(std::ostream& out)
   { }
  
   /*
   * Clear all timers.
   */ 
   template<int D>
   void Simulator<D>::clearTimers()
   { }

}
}
#endif
