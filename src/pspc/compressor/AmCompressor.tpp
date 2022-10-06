#ifndef PSPC_AM_COMPRESSOR_TPP
#define PSPC_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <pspc/System.h>
#include <util/global.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   // Constructor
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
   : Compressor<D>(system)
   {  setClassName("AmCompressor"); }

   // Destructor
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {  setClassName("AmCompressor"); }

   // Read parameters from file
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>, DArray<double> >::readParameters(in);
   }

   // Compute and return L2 norm of residual vector
   template <int D>
   double AmCompressor<D>::findNorm(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double normResSq = 0.0;
      for (int i = 0; i < n; i++) {
         normResSq += hist[i] * hist[i];
      }
      return sqrt(normResSq);
   }

   // Compute and return maximum element of residual vector.
   template <int D>
   double AmCompressor<D>::findMaxAbs(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double maxRes = 0.0;
      for (int i = 0; i < n; i++) {
         if (fabs(hist[i]) > maxRes)
            maxRes = fabs(hist[i]);
      }
      return maxRes;
   }

   // Update basis
   template <int D>
   void AmCompressor<D>::updateBasis(RingBuffer<DArray<double> > & basis,
                                   RingBuffer<DArray<double> > const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      DArray<double> newbasis;
      newbasis.allocate(n);

      for (int i = 0; i < n; i++) {
         // sequential histories basis vectors
         newbasis[i] = hists[0][i] - hists[1][i]; 
      }

      basis.append(newbasis);
   }

   // Compute one element of U matrix of by computing a dot product
   template <int D>
   double
   AmCompressor<D>::computeUDotProd(RingBuffer<DArray<double> > const & resBasis,
                                  int m, int n)
   {
      const int length = resBasis[0].capacity();
      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resBasis[m][i] * resBasis[n][i];
      }
      return dotprod;
   }

   // Compute one element of V vector by computing a dot product
   template <int D>
   double
   AmCompressor<D>::computeVDotProd(DArray<double> const & resCurrent,
                                  RingBuffer<DArray<double> > const & resBasis,
                                  int m)
   {
      const int length = resBasis[0].capacity();
      double dotprod = 0.0;
      for(int i = 0; i < length; i++) {
         dotprod += resCurrent[i] * resBasis[m][i];
      }
      return dotprod;
   }

   // Update entire U matrix
   template <int D>
   void AmCompressor<D>::updateU(DMatrix<double> & U,
                               RingBuffer<DArray<double> > const & resBasis,
                               int nHist)
   {
      // Update matrix U by shifting elements diagonally
      int maxHist = U.capacity1();
      for (int m = maxHist-1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1);
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist; ++m) {
         double dotprod = computeUDotProd(resBasis,0,m);
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   template <int D>
   void AmCompressor<D>::updateV(DArray<double> & v,
                               DArray<double> const & resCurrent,
                               RingBuffer<DArray<double> > const & resBasis,
                               int nHist)
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }

   template <int D>
   void AmCompressor<D>::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

   template <int D>
   void
   AmCompressor<D>::addHistories(DArray<double>& trial,
                               RingBuffer<DArray<double> > const & basis,
                               DArray<double> coeffs,
                               int nHist)
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist; i++) {
         for (int j = 0; j < n; j++) {
            // Not clear on the origin of the -1 factor
            trial[j] += coeffs[i] * -1 * basis[i][j];
         }
      }
   }

   template <int D>
   void AmCompressor<D>::addPredictedError(DArray<double>& fieldTrial,
                                         DArray<double> const & resTrial,
                                         double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }

   // Does the system have an initial field guess?
   template <int D>
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   // Compute and return the number of elements in a field vector
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().basis().nBasis(); }

   // Get the current field from the system
   template <int D>
   void AmCompressor<D>::getCurrent(DArray<double>& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      // Application of the template is not straightforward.

   }

   // Perform the main system computation (solve the MDE)
   template <int D>
   void AmCompressor<D>::evaluate()
   {  system().compute(); }

   // Compute the residual for the current system state
   template <int D>
   void AmCompressor<D>::getResidual(DArray<double>& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      // Initialize residuals
      for (int i = 1 ; i < n; ++i) {
         resid[i] = 0.0;
      }
      resid[0] = -1.0;

       // Compute SCF residual vector elements
      for (int j = 0; j < nMonomer; ++j) {
        for (int k = 0; k < nBasis; ++k) {
           resid[k] += system().c().basis(j)[k];
        }
      }

   }

   // Update the current system field coordinates
   template <int D>
   void AmCompressor<D>::update(DArray<double>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      // Need a strategy for this - it's not straightforward

   }

   template<int D>
   void AmCompressor<D>::outputToLog()
   {}

}
}
#endif
