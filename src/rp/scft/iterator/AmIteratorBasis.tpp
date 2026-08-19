#ifndef RP_AM_ITERATOR_BASIS_TPP
#define RP_AM_ITERATOR_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorBasis.h"

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>

#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/containers/DArray.h>

#include <util/global.h>
#include <cmath>

#include <pscf/iterator/AmIteratorTmpl.tpp> // base template implementation

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D, class T>
   AmIteratorBasis<D,T>::AmIteratorBasis(System<D,T>& system)
    : AmIterTmplT(),
      interaction_(),
      scaleStress_(10.0)
   {
      Iterator<D,T>::setSystem(system);
      Iterator<D,T>::isSymmetric_ = true;  
      ParamComposite::setClassName("AmIteratorBasis"); 
   }

   /*
   * Destructor.
   */
   template <int D, class T>
   AmIteratorBasis<D,T>::~AmIteratorBasis()
   {}

   /*
   * Read parameters from file.
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::readParameters(std::istream& in)
   {
      // Use base class methods to read parameters
      AmIterTmplT::readParameters(in);
      AmIterTmplT::readErrorType(in);

      UnitCell<D> const & unitCell = system().domain().unitCell();
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      int np = unitCell.nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);

      // Read optional isFlexible boolean (true by default)
      Iterator<D,T>::isFlexible_ = 1;  // Default
      ParamComposite::readOptional(in, "isFlexible", 
                                   Iterator<D,T>::isFlexible_);

      // Populate flexibleParams_ bool array
      if (Iterator<D,T>::isFlexible_) {
         flexibleParams_.clear();
         // Set all flexibleParams_ values to true by default
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(true); 
         }
         // Optionally read flexibleParams_ array
         ParamComposite::readOptionalFSArray(in, "flexibleParams", 
                                             flexibleParams_, np);
         if (Iterator<D,T>::nFlexibleParams() == 0) {
            Iterator<D,T>::isFlexible_ = false;
         }
      } else { // If isFlexible_ == false
         // Set all flexibleParams_ values to false
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(false);
         }
      }

      // Optionally read scaleStress value
      scaleStress_ = 10.0; // Default value
      ParamComposite::readOptional(in, "scaleStress", scaleStress_);

      // Read optional mixing parameters (lambda, useLambdaRamp, and r)
      AmIterTmplT::readMixingParameters(in);

      // Allocate local AmbdInteraction class
      const int nMonomer = system().mixture().nMonomer();
      interaction_.setNMonomer(nMonomer);

   }

   /*
   * Output timing results to log file.
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Iterator times contributions:\n";
      AmIterTmplT::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D, class T>
   void AmIteratorBasis<D,T>::setup(bool isContinuation)
   {
      AmIterTmplT::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions

   /*
   * Compute and return number of elements in a residual vector.
   */
   template <int D, class T>
   int AmIteratorBasis<D,T>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      int nEle = nMonomer*nBasis;
      if (Iterator<D,T>::isFlexible_) {
         nEle += Iterator<D,T>::nFlexibleParams();
      }
      return nEle;
   }

   /*
   * Does the system have an initial field guess?
   */
   template <int D, class T>
   bool AmIteratorBasis<D,T>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current state (w fields and lattice parameters).
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::getCurrent(VectorT& state)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const int nEle = nElements();
      UTIL_CHECK(state.capacity() == nEle);

      // Copy all w-fields into a linear array
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         DArray<double> const & field = system().w().basis(i);
         begin = i*nBasis;
         for (int k = 0; k < nBasis; k++) {
            state[begin + k] = field[k];
         }
      }

      // Add elements associated with unit cell parameters (if any)
      if (Iterator<D,T>::isFlexible_) {
         UTIL_CHECK((Iterator<D,T>::nFlexibleParams() > 0));
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double,6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         begin = nMonomer*nBasis;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               state[begin + counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK((counter == Iterator<D,T>::nFlexibleParams()));
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::evaluate()
   {  system().compute(Iterator<D,T>::isFlexible_); }

   /*
   * Compute the residual for the current system state.
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::getResidual(VectorT& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();
      const int n = nElements();
      UTIL_CHECK(resid.capacity() == n);

      // Local variables
      double chi, p;
      int i, j, k, begin;

      // Initialize residual vector to zero
      for (i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (i = 0; i < nMonomer; ++i) {
         begin = i*nBasis;
         for (j = 0; j < nMonomer; ++j) {
            chi = interaction_.chi(i,j);
            p = interaction_.p(i,j);
            DArray<double> const & c = system().c().basis(j);
            DArray<double> const & w = system().w().basis(j);
            if (system().h().hasData()) {
               DArray<double> const & h = system().h().basis(j);
               for (k = 0; k < nBasis; ++k) {
                  resid[begin + k] += chi*c[k] + p*(h[k] - w[k]);
               }
            } else {
               for (k = 0; k < nBasis; ++k) {
                  resid[begin + k] += chi*c[k] - p*w[k];
               }
            }
         }
      }

      // Add term proportional to sum of all concentrations
      double shift = -1.0 / interaction_.sumChiInverse();
      if (system().mask().hasData()) {
         DArray<double> const & mask = system().mask().basis();
         for (i = 0; i < nMonomer; ++i) {
            begin = i*nBasis;
            for (k = 0; k < nBasis; ++k) {
               resid[begin + k] += shift*mask[k];
            }
         }
      } else {
         for (i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] += shift;
         }
      }

      // If canonical ensemble, zero homogeneous residual components
      if (system().mixture().isCanonical()) {
         for (i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If flexible unit cell, then compute stress residuals
      if (Iterator<D,T>::isFlexible_) {

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         double coeff = -1.0 * scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         begin = nMonomer*nBasis;
         int counter = 0;
         for (i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               resid[begin + counter] = coeff * Iterator<D,T>::stress(i);
               counter++;
            }
         }
         UTIL_CHECK((counter == Iterator<D,T>::nFlexibleParams()));
      }

   }

   /*
   * Update the current system field vector.
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::update(VectorT& newState)
   {
      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().domain().basis().nBasis();

      // Allocate wFields container
      DArray< DArray<double> > wFields;
      wFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         wFields[i].allocate(nBasis);
      }

      // Copy w fields from newState to wFields container
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         begin = i*nBasis;
         for (int k = 0; k < nBasis; k++) {
            wFields[i][k] = newState[begin + k];
         }
      }

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {

         // Set homogeneous components of all w fields to zero
         for (int i = 0; i < nMonomer; ++i) {
            wFields[i][0] = 0.0;
         }

         // Add average values arising from interactions
         double chi, wAve, cAve;
         for (int i = 0; i < nMonomer; ++i) {
            wAve = 0.0;
            for (int j = 0; j < nMonomer; ++j) {
               chi = interaction_.chi(i,j);
               cAve = system().c().basis(j)[0];
               wAve += chi * cAve;
            }
            wFields[i][0] = wAve;
         }

         // If external fields exist, add their spatial averages
         if (system().h().hasData()) {
            for (int i = 0; i < nMonomer; ++i) {
               wFields[i][0] += system().h().basis(i)[0];
            }
         }
      }

      // Set fields in system w container
      system().w().setBasis(wFields);

      // If flexible, update unit cell parameters
      if (Iterator<D,T>::isFlexible_) {

         // Initialize parameters array with current values
         FSArray<double, 6> parameters;
         parameters = system().domain().unitCell().parameters();

         // Reset any parameters that are flexible
         const int nParam = system().domain().unitCell().nParameter();
         const int begin = nMonomer*nBasis;
         const double scale = 1.0 / scaleStress_;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = scale * newState[begin + counter];
               counter++;
            }
         }
         UTIL_CHECK((counter == Iterator<D,T>::nFlexibleParams()));

         // Set system unit cell parameters
         system().setUnitCell(parameters);
      }
   }

   /*
   * Output relevant system details to the iteration log file.
   */
   template <int D, class T>
   void AmIteratorBasis<D,T>::outputToLog()
   {
      if (Iterator<D,T>::isFlexible_ && AmIterTmplT::verbose() > 1) {

         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nBasis = system().domain().basis().nBasis();
         const int begin = nMonomer*nBasis;

         double res, str;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               res = AmIterTmplT::residual()[begin + counter];
               str = -1.0 * res / scaleStress_;
               Log::file() 
                  << " Cell Param  " << i << " = "
                  << Dbl(unitCell.parameters()[i], 15)
                  << " , stress = " 
                  << Dbl(str, 15)
                  << "\n";
               counter++;
            }
         }
      }
   }

}
}
#endif
