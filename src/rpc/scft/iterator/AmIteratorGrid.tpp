#ifndef RPC_AM_ITERATOR_GRID_TPP
#define RPC_AM_ITERATOR_GRID_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/cpu/Reduce.h>
#include <prdc/cpu/VecOp.h>
#include <pscf/inter/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/containers/RingBuffer.h>
#include <util/containers/DArray.h>
#include <util/global.h>
#include <cmath>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
    : Iterator<D>(system)
   {
      isSymmetric_ = true;
      ParamComposite::setClassName("AmIteratorGrid");
   }

   /*
   * Destructor.
   */
   template <int D>
   AmIteratorGrid<D>::~AmIteratorGrid()
   {}

   /*
   * Read parameters from file.
   */
   template <int D>
   void AmIteratorGrid<D>::readParameters(std::istream& in)
   {
      // Use base class methods to read parameters
      AmTmpl::readParameters(in);
      AmTmpl::readErrorType(in);

      UnitCell<D> const & unitCell = system().domain().unitCell();
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      int np = unitCell.nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);

      // Read optional isFlexible boolean (true by default)
      isFlexible_ = 1;  // Default
      readOptional(in, "isFlexible", isFlexible_);

      // Populate flexibleParams_ based on isFlexible_ (all 0s or all 1s),
      // then optionally overwrite with user input from param file
      if (isFlexible_) {
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(true); // Set all values to true
         }
         // Read optional flexibleParams_ array to overwrite current array
         readOptionalFSArray(in, "flexibleParams", flexibleParams_, np);
         if (nFlexibleParams() == 0) {
            isFlexible_ = false;
         }
      } else { // isFlexible_ = false
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(false); // Set all values to false
         }
      }

      // Read optional scaleStress value
      scaleStress_ = 10.0; // Default
      readOptional(in, "scaleStress", scaleStress_);

      // Read option mixing parameters (lambda, useLambdaRamp, and r)
      AmTmpl::readMixingParameters(in);

      // Allocate local modified copy of Interaction class
      const int nMonomer = system().mixture().nMonomer();
      interaction_.setNMonomer(nMonomer);

   }

   /*
   * Output timing results to log file.
   */
   template<int D>
   void AmIteratorGrid<D>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Iterator times contributions:\n";
      AmTmpl::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorGrid<D>::setup(bool isContinuation)
   {
      AmTmpl::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions to exchange data with parent system

   /*
   * Compute the number of elements in the residual vector.
   */
   template <int D>
   int AmIteratorGrid<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();

      int nEle = nMonomer*nMesh;
      if (isFlexible_) {
         nEle += nFlexibleParams();
      }
      return nEle;
   }

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool AmIteratorGrid<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field vector (w fields and lattice parameters).
   */
   template <int D>
   void AmIteratorGrid<D>::getCurrent(DArray<double>& state)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const int nEle = nElements();
      UTIL_CHECK(state.capacity() == nEle);

      // Copy w-fields into a linear array
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         const RField<D>& field = system().w().rgrid(i);
         begin = i*nMesh;
         for (int k = 0; k < nMesh; k++) {
            state[begin + k] = field[k];
         }
      }

      // Add elements associated with unit cell parameters (if any)
      if (isFlexible()) {
         UTIL_CHECK(nFlexibleParams() > 0);
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double, 6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         begin = nMonomer*nMesh;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               state[begin + counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void AmIteratorGrid<D>::evaluate()
   {  system().compute(isFlexible_); }

   /*
   * Compute the residual for the current system state.
   */
   template <int D>
   void AmIteratorGrid<D>::getResidual(DArray<double>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
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
         begin = i*nMesh;
         for (j = 0; j < nMonomer; ++j) {
            chi = interaction_.chi(i,j);
            p = interaction_.p(i,j);
            RField<D> const & c = system().c().rgrid(j);
            RField<D> const & w = system().w().rgrid(j);
            if (system().h().hasData()) {
               RField<D> const & h = system().h().rgrid(j);
               for (k = 0; k < nMesh; ++k) {
                  resid[begin + k] += chi*c[k] + p*(h[k] - w[k]) ;
               }
            } else {
               for (k = 0; k < nMesh; ++k) {
                  resid[begin + k] += chi*c[k] - p*w[k];
               }
            }
         }
      }

      // Add term proportional to sum of all concentrations
      double shift = -1.0 / interaction_.sumChiInverse();
      if (system().mask().hasData()) {
         RField<D> const & mask = system().mask().rgrid();
         for (i = 0; i < nMonomer; ++i) {
            begin = i*nMesh;
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] += shift*mask[k];
            }
         }
      } else {
         for (i = 0; i < nMonomer; ++i) {
            begin = i*nMesh;
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] += shift;
            }
         }
      }

      // If system is canonical, set all spatial averages to zero
      if (system().mixture().isCanonical()) {
         double average;
         for (i = 0; i < nMonomer; ++i) {
            begin = i*nMesh;
            average = 0.0;
            for (k = 0; k < nMesh; ++k) {
               average += resid[begin + k];
            }
            average /= double(nMesh);
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] -= average;
            }
         }
      }

      // If flexible unit cell, then compute stress residuals
      if (isFlexible()) {

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         const double coeff = -1.0 * scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         begin = nMonomer*nMesh;
         int counter = 0;
         for (i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               resid[begin + counter] = coeff * stress(i);
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());
      }

   }

   /*
   * Update the current system field vector.
   */
   template <int D>
   void AmIteratorGrid<D>::update(DArray<double>& newState)
   {

      // References to system components
      Domain<D> const & domain = system().domain();
      Mesh<D> const & mesh = domain.mesh();
      IntVec<D> const & meshDimensions = mesh.dimensions();

      // Constants
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = mesh.size();

      // Allocate wFields container
      DArray< RField<D> > wFields;
      wFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         wFields[i].allocate(meshDimensions);
         UTIL_CHECK(wFields[i].capacity() == nMesh);
      }

      // Copy new fields from newState vector to wFields container
      int begin;
      for (int i = 0; i < nMonomer; i++) {
         begin = i*nMesh;
         for (int k = 0; k < nMesh; k++) {
            wFields[i][k] = newState[begin + k];
         }
      }

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {

         // Subtract spatial average from each w field
         double wAve;
         for (int i = 0; i < nMonomer; ++i) {
            wAve = Reduce::sum(wFields[i]);
            wAve /= double(nMesh);
            VecOp::subEqS(wFields[i], wAve);
         }

         // Compute spatial averages of all concentration fields
         DArray<double> cAve;
         cAve.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            cAve[i] = Reduce::sum(system().c().rgrid(i));
            cAve[i] /= double(nMesh);
         }

         // Add average value arising from interactions
         double chi;
         for (int i = 0; i < nMonomer; ++i) {
            wAve = 0.0;
            for (int j = 0; j < nMonomer; ++j) {
               chi = interaction_.chi(i,j);
               wAve += chi * cAve[j];
            }
            VecOp::addEqS(wFields[i], wAve);
         }

         // If external fields exist, add their spatial averages
         if (system().h().hasData()) {
            double hAve;
            for (int i = 0; i < nMonomer; ++i) {
               hAve = Reduce::sum(system().h().rgrid(i));
               hAve /= double(nMesh);
               VecOp::addEqS(wFields[i], hAve);
            }
         }
      }

      // Set fields in system w container
      system().w().setRGrid(wFields);

      // If flexible, update unit cell parameters
      if (isFlexible()) {

         // Initialize parameters array with current values
         FSArray<double, 6> parameters;
         parameters = domain.unitCell().parameters();

         // Reset any parameters that are flexible
         const int nParam = domain.unitCell().nParameter();
         double coeff = 1.0 / scaleStress_;
         const int begin = nMonomer*nMesh;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = coeff * newState[begin + counter];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlexibleParams());

         // Set system unit cell parameters
         system().setUnitCell(parameters);
      }

   }

   /*
   * Output relevant system details to the iteration log.
   */
   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible() && verbose() > 1) {
         double str;
         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nMesh = system().domain().mesh().size();
         const int begin = nMonomer*nMesh;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               str = - 1.0 * residual()[begin + counter] / scaleStress_;
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
