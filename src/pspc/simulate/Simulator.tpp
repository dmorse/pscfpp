#ifndef PSPC_SIMULATOR_TPP
#define PSPC_SIMULATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulator.h"
#include <pspc/System.h>
#include <pspc/compressor/Compressor.h>
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
   Simulator<D>::Simulator(System<D>& system)
   : random_(),
     iStep_(0),
     hasHamiltonian_(false),
     hasWC_(false),
     systemPtr_(&system)
   { setClassName("Simulator"); }

   /*
   * Destructor.
   */
   template <int D>
   Simulator<D>::~Simulator()
   {}

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void Simulator<D>::readParameters(std::istream &in)
   {
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
   * Perform a field theoretic MC simulation of nStep steps.
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
   * Compute Monte Carlo Hamiltonian.
   */
   template <int D>
   void Simulator<D>::computeHamiltonian()
   {
      UTIL_CHECK(system().w().hasData());
      UTIL_CHECK(system().hasCFields());
      UTIL_CHECK(hasWC_);
      hasHamiltonian_ = false;

      Mixture<D> const & mixture = system().mixture();
      Domain<D> const & domain = system().domain();

      const int nMonomer = mixture.nMonomer();
      const int meshSize = domain.mesh().size();

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
      // lnQ now contains a value per monomer

      // Compute field contribution HW
      int i, j;
      double prefactor, w;
      double HW = 0.0;
      // Compute quadratic field contribution to HW
      for (j = 0; j < nMonomer - 1; ++j) {
         prefactor = -0.5*double(nMonomer)/chiEvals_[j];
         RField<D> const & wc = wc_[j];
         for (i = 0; i < meshSize; ++i) {
            w = wc[i];
            HW += prefactor*w*w;
         }
      }
      
      // Subtract average of Langrange multiplier field
      RField<D> const & xi = wc_[nMonomer-1];
      for (i = 0; i < meshSize; ++i) {
         HW -= xi[i];
      }

      // Normalize HW to equal a value per monomer
      HW /= double(meshSize);
      
      // Compute final MC Hamiltonian
      const double vSystem  = domain.unitCell().volume();
      const double vMonomer = mixture.vMonomer();
      const double nMonomerSystem = vSystem / vMonomer;
      fieldHamiltonian_ = nMonomerSystem * HW;
      idealHamiltonian_ = -1.0 * nMonomerSystem * lnQ;
      hamiltonian_ = idealHamiltonian_ + fieldHamiltonian_;

      hasHamiltonian_ = true;
   }

   template <int D>
   void Simulator<D>::analyzeChi()
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
   void Simulator<D>::computeWC()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;
      // Loop over eigenvectors (j is an eigenvector index)
      for (j = 0; j < nMonomer; ++j) {
         // Loop over grid points to zero out field wc_[j]
         RField<D>& wc = wc_[j];
         for (i = 0; i < meshSize; ++i) {
            wc[i] = 0.0;
         }

         // Loop over monomer types (k is a monomer index)
         for (k = 0; k < nMonomer; ++k) {
            double vec = chiEvecs_(j, k)/double(nMonomer);

            // Loop over grid points
            RField<D> const & w = system().w().rgrid(k);
            for (i = 0; i < meshSize; ++i) {
               wc[i] += vec*w[i];
            }

         }
      }
      
      // Debugging output
      #if 0
      Log::file() << "wc " << wc_.capacity() << "\n";
      for (i = 0; i < 10; ++i) {
         Log::file() << "wc_1 " << wc_[0][i] << "\n";
         Log::file() << "wc_2 " << wc_[1][i] << "\n";
      }
      #endif

      hasWC_ = true;
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