/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIterator.h"
#include <r1d/System.h>
#include <pscf/inter/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/global.h>

namespace Pscf {
namespace R1d{

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   AmIterator::AmIterator(System& system)
   : Iterator(system),
     interaction_()
   {  setClassName("AmIterator"); }

   /*
   * Destructor.
   */
   AmIterator::~AmIterator()
   {}

   /*
   * Read parameters from file.
   */
   void AmIterator::readParameters(std::istream& in)
   {
      // Call parent class readParameters and readErrorType functions
      AmIteratorTmpl<Iterator, DArray<double> >::readParameters(in);
      AmIteratorTmpl<Iterator, DArray<double> >::readErrorType(in);

      const int nMonomer = system().mixture().nMonomer(); 
      interaction_.setNMonomer(nMonomer);
   }

   // Protected member function

   /*
   * Setup before entering iteration loop.
   */
   void AmIterator::setup(bool isContinuation)
   {
      AmIteratorTmpl<Iterator, DArray<double> >::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual member functions

   // Compute and return the number of elements in a field vector
   int AmIterator::nElements()
   {
      const int nm = system().mixture().nMonomer(); // # of monomers
      const int nx = domain().nx();                 // # of grid points
      return nm*nx;
   }

   bool AmIterator::hasInitialGuess()
   {
      // R1d::System doesn't have an appropriate query function
      return true;
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

   /*
   * Solve the MDEs for the current system w fields.
   */
   void AmIterator::evaluate()
   {  mixture().compute(system().wFields(), system().cFields()); }

   /*
   * Check whether this is canonical ensemble.
   */
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
      const double shift = -1.0/interaction_.sumChiInverse();
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
            chi = interaction_.chi(i,j);
            p   = interaction_.p(i,j);
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
