/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <fd1d/System.h>
#include <pscf/inter/Interaction.h>
#include <util/global.h>

namespace Pscf {
namespace Fd1d{

   using namespace Util;

   // Constructor
   AmIterator::AmIterator(System& system)
   : Iterator(system)
   {  setClassName("AmIterator"); }

   // Destructor
   AmIterator::~AmIterator()
   {}

   // Read parameters from file
   void AmIterator::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Iterator, DArray<double> >::readParameters(in);
   }

   // Compute and return L2 norm of residual vector
   double AmIterator::findNorm(DArray<double> const & hist)
   {
      const int n = hist.capacity();
      double normResSq = 0.0;
      for (int i = 0; i < n; i++) {
         normResSq += hist[i] * hist[i];
      }
      return sqrt(normResSq);
   }

   // Compute and return maximum element of residual vector.
   double AmIterator::findMaxAbs(DArray<double> const & hist)
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
   void AmIterator::updateBasis(RingBuffer<DArray<double> > & basis,
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
   double
   AmIterator::computeUDotProd(RingBuffer< DArray<double> > const& resBasis,
                               int m, int n)
   {
      const int nr = resBasis[0].capacity();
      double dotprod = 0.0;
      for (int i = 0; i < nr; i++) {
         dotprod += resBasis[m][i] * resBasis[n][i];
      }
      return dotprod;
   }

   // Compute one element of V vector by computing a dot product
   double
   AmIterator::computeVDotProd(DArray<double> const & resCurrent,
                               RingBuffer<DArray<double> > const& resBasis,
                               int m)
   {
      const int nr = resBasis[0].capacity();
      double dotprod = 0.0;
      for (int i = 0; i < nr; i++) {
         dotprod += resCurrent[i] * resBasis[m][i];
      }
      return dotprod;
   }

   // Update entire U matrix
   void AmIterator::updateU(DMatrix<double> & U,
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

   void AmIterator::updateV(DArray<double> & v,
                            DArray<double> const & resCurrent,
                            RingBuffer<DArray<double> > const & resBasis,
                            int nHist)
   {
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }

   void AmIterator::setEqual(DArray<double>& a, DArray<double> const & b)
   {  a = b; }

   void AmIterator::addHistories(DArray<double>& trial,
                                 RingBuffer< DArray<double> > const& basis,
                                 DArray<double> coeffs,
                                 int nHist)
   {
      int n = trial.capacity();
      for (int i = 0; i < nHist; i++) {
         for (int j = 0; j < n; j++) {
            trial[j] += coeffs[i] * -1 * basis[i][j];
         }
      }
   }

   void AmIterator::addPredictedError(DArray<double>& fieldTrial,
                                      DArray<double> const & resTrial,
                                      double lambda)
   {
      int n = fieldTrial.capacity();
      for (int i = 0; i < n; i++) {
         fieldTrial[i] += lambda * resTrial[i];
      }
   }

   bool AmIterator::hasInitialGuess()
   {
      // Fd1d::System doesn't hav a hasFields() function
      return true;
   }

   // Compute and return the number of elements in a field vector
   int AmIterator::nElements()
   {
      const int nm = system().mixture().nMonomer(); // # of monomers
      const int nx = domain().nx();                 // # of grid points
      return nm*nx;
   }

   // Get the current fields from the system as a 1D array
   void AmIterator::getCurrent(DArray<double>& curr)
   {
      const int nm = system().mixture().nMonomer();  // # of monomers
      const int nx =  domain().nx();                 // # of grid points
      const DArray< DArray<double> > * currSys = &system().wFields();

      // Straighten out fields into linear arrays
      for (int i = 0; i < nm; i++) {
         for (int k = 0; k < nx; k++) {
            curr[i*nx+k] = (*currSys)[i][k];
         }
      }
   }

   // Solve the MDEs for the current system w fields
   void AmIterator::evaluate()
   {
      mixture().compute(system().wFields(), system().cFields());
   }


   // Check whether Canonical ensemble
   bool AmIterator::isCanonical()
   {
      Species::Ensemble ensemble;

      // Check ensemble of all polymers
      for (int i = 0; i < mixture().nPolymer(); ++i) {
         ensemble = mixture().polymer(i).ensemble();
         if (ensemble == Species::Open) {
            return false;
         }
         if (ensemble == Species::Unknown) {
            UTIL_THROW("Unknown species ensemble");
         }
      }

      // Check ensemble of all solvents
      for (int i = 0; i < mixture().nSolvent(); ++i) {
         ensemble = mixture().solvent(i).ensemble();
         if (ensemble == Species::Open) {
            return false;
         }
         if (ensemble == Species::Unknown) {
            UTIL_THROW("Unknown species ensemble");
         }
      }

      // Return true if no species had an open ensemble
      return true;
   }

   // Compute the residual for the current system state
   void AmIterator::getResidual(DArray<double>& resid)
   {
      const int nm = system().mixture().nMonomer();
      const int nx = domain().nx();
      const int nr = nm*nx;

       // Initialize residuals
      const double shift = -1.0/system().interaction().sum_inv();
      for (int i = 0 ; i < nr; ++i) {
         resid[i] = shift;
      }

      // Compute SCF residual vector elements
      double chi, p;
      int i, j, k;
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nm; ++j) {
            DArray<double>& cField = system().cField(j);
            DArray<double>& wField = system().wField(j);
            chi = system().interaction().chi(i,j);
            p   = system().interaction().idemp(i,j);
            for (k = 0; k < nx; ++k) {
               int idx = i*nx + k;
               resid[idx] += chi*cField[k] - p*wField[k];
            }
         }
      }

   }

   // Update the current system field coordinates
   void AmIterator::update(DArray<double>& newGuess)
   {
      const int nm = mixture().nMonomer();  // # of monomers
      const int nx = domain().nx();         // # of grid points

      // Set current w fields in system
      // If canonical, shift to as to set last element to zero
      const double shift = newGuess[nm*nx - 1];
      for (int i = 0; i < nm; i++) {
         DArray<double>& wField = system().wField(i);
         if (isCanonical()) {
            for (int k = 0; k < nx; k++) {
               wField[k] = newGuess[i*nx + k] - shift;
            }
         } else {
            for (int k = 0; k < nx; k++) {
               wField[k] = newGuess[i*nx + k];
            }
         }
      }

   }

   void AmIterator::outputToLog()
   {}

}
}
