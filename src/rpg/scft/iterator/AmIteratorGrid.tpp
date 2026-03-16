#ifndef RPG_AM_ITERATOR_GRID_TPP
#define RPG_AM_ITERATOR_GRID_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <prdc/cuda/resources.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/NanException.h>
#include <util/global.h>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
    : Iterator<D>(system)
   {
      ParamComposite::setClassName("AmIteratorGrid");
      Iterator<D>::isSymmetric_ = false;
   }

   /*
   * Destructor.
   */
   template <int D>
   AmIteratorGrid<D>::~AmIteratorGrid()
   {}

   /*
   * Read parameter file block.
   */
   template <int D>
   void AmIteratorGrid<D>::readParameters(std::istream& in)
   {
      // Preconditions on unit cell
      UnitCell<D> const & unitCell = system().domain().unitCell();
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      int np = unitCell.nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);

      // Read param file format for base class
      AmIterTmplT::readParameters(in);
      AmIterTmplT::readErrorType(in);

      // Read optional isFlexible boolean (true by default)
      isFlexible_ = 1;
      ParamComposite::readOptional(in, "isFlexible", isFlexible_);

      // Populate flexibleParams_ based on isFlexible_ (all 0s or all 1s),
      // then optionally overwrite with user input from param file
      if (isFlexible_) {
         flexibleParams_.clear();
         for (int i = 0; i < np; ++i) {
            flexibleParams_.append(true); // Set all values to true
         }
         // Read optional flexibleParams_ array to overwrite current array
         ParamComposite::readOptionalFSArray(in, "flexibleParams",
                                             flexibleParams_, np);
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
      scaleStress_ = 10.0;  // default
      ParamComposite::readOptional(in, "scaleStress", scaleStress_);

      // Read optional mixing parameters (lambda, useLambdaRamp, r)
      AmIterTmplT::readMixingParameters(in);

      // Allocate local modified copy of Interaction class
      interaction_.setNMonomer(system().mixture().nMonomer());
   }

   /*
   * Output timing results to log file.
   */
   template<int D>
   void AmIteratorGrid<D>::outputTimers(std::ostream& out) const
   {
      out << "\n";
      out << "Iterator times contributions:\n";
      AmIterTmplT::outputTimers(out);
   }

   // Protected virtual function

   /*
   * Setup before entering iteration loop.
   */
   template <int D>
   void AmIteratorGrid<D>::setup(bool isContinuation)
   {
      AmIterTmplT::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions

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
   * Check if the system has an initial guess.
   */
   template <int D>
   bool AmIteratorGrid<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current w fields and lattice parameters.
   */
   template <int D>
   void AmIteratorGrid<D>::getCurrent(VectorT& state)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const int n = nElements();
      UTIL_CHECK(state.capacity() == n);

      // Copy all system w-fields into a linear array
      VectorT slice;
      for (int i = 0; i < nMonomer; i++) {
         slice.associate(state, i*nMesh, nMesh);
         VecOp::eqV(slice, system().w().rgrid(i));
         slice.dissociate();
      }

      // If flexible unit cell, also store unit cell parameters
      if (isFlexible_) {
         UTIL_CHECK(nFlexibleParams() > 0);
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double, 6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         DArray<cudaReal> paramsTmp(nFlexibleParams());
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               paramsTmp[counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == paramsTmp.capacity());

         // Copy unit cell parameters to the end of the state array
         slice.associate(state, nMonomer*nMesh, paramsTmp.capacity());
         slice = paramsTmp; // copy from host to device, for GPU code
         slice.dissociate();
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
   void AmIteratorGrid<D>::getResidual(VectorT& resid)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const int n = nElements();
      UTIL_CHECK(resid.capacity() == n);

      // Initialize residual vector to zero.
      VecOp::eqS(resid, 0.0);

      // Array of VectorT arrays associated with slices of resid.
      // one VectorT array per monomer species, each of size nMesh.
      DArray<VectorT> residSlices;
      residSlices.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         residSlices[i].associate(resid, i*nMesh, nMesh);
      }

      // Compute SCF residuals
      for (int i = 0; i < nMonomer; i++) {
         for (int j = 0; j < nMonomer; j++) {
            VecOp::addVcVcVc(residSlices[i], residSlices[i], 1.0,
                             system().c().rgrid(j), interaction_.chi(i, j),
                             system().w().rgrid(j), -interaction_.p(i, j));
         }
      }

      // If iterator has mask, account for it in residual values
      if (system().mask().hasData()) {
         double coeff = -1.0 / interaction_.sumChiInverse();
         for (int i = 0; i < nMonomer; ++i) {
            VecOp::addEqVc(residSlices[i], system().mask().rgrid(), 
                           coeff);
         }
      }

      // If iterator has external fields, account for them in the values
      // of the residuals
      if (system().h().hasData()) {
         for (int i = 0; i < nMonomer; ++i) {
            for (int j = 0; j < nMonomer; ++j) {
               double p = interaction_.p(i,j);
               VecOp::addEqVc(residSlices[i], system().h().rgrid(j), p);
            }
         }
      }

      // If ensemble is not canonical, account for incompressibility.
      if (!system().mixture().isCanonical()) {
         cudaReal factor = 1.0 / interaction_.sumChiInverse();
         VecOp::subEqS(resid, factor, 0, nMonomer*nMesh);
      } else {
         for (int i = 0; i < nMonomer; i++) {
            // Find current average
            cudaReal average = findAverage(residSlices[i]);
            // subtract out average to set residual average to zero
            VecOp::subEqS(residSlices[i], average);
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible_) {

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
         HostDArray<cudaReal> stressH(nFlexibleParams());
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               stressH[counter] = coeff * Iterator<D>::stress(i);
               counter++;
            }
         }
         UTIL_CHECK(counter == stressH.capacity());

         VectorT stressD;
         stressD.associate(resid, nMonomer*nMesh, stressH.capacity());
         stressD = stressH; // copy from host to device
         UTIL_CHECK(stressD.isAssociated());
         stressD.dissociate();
      }

      // Dissociate elements of residSlices from resid array
      for (int i = 0; i < nMonomer; i++) {
         UTIL_CHECK(residSlices[i].isAssociated());
         residSlices[i].dissociate();
      }

   }

   /*
   * Update the system with a new trial field vector.
   */
   template <int D>
   void AmIteratorGrid<D>::update(VectorT& newGuess)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {
         cudaReal average, wAverage, cAverage;
         for (int i = 0; i < nMonomer; i++) {

            // Define array associated with a slice of newGuess
            VectorT ngSlice;
            ngSlice.associate(newGuess, i*nMesh, nMesh);

            // Find current spatial average
            average = findAverage(ngSlice);

            // Subtract average from field, setting average to zero
            VecOp::subEqS(ngSlice, average);

            // Compute the average omega(i) due to interactions
            wAverage = 0.0;
            for (int j = 0; j < nMonomer; j++) {
               cAverage = findAverage(system().c().rgrid(j));
               wAverage += interaction_.chi(i,j) * cAverage;
            }

            // If external fields exist, add them to wAverage
            if (system().h().hasData()) {
               wAverage += findAverage(system().h().rgrid(i));
            }

            // Add new average omega value to the field
            VecOp::addEqS(ngSlice, wAverage);
            ngSlice.dissociate();
         }
      }

      // Set w fields in system w container
      system().w().setRGrid(newGuess);

      // If flexible unit cell, update cell parameters
      if (isFlexible_) {
         FSArray<double,6> parameters;
         parameters = system().domain().unitCell().parameters();
         const int nParam = system().domain().unitCell().nParameter();

         // transfer from device to host
         HostDArray<cudaReal> tempH(nFlexibleParams());
         tempH.copySlice(newGuess, nMonomer*nMesh);

         const double coeff = 1.0 / scaleStress_;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = coeff * tempH[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == tempH.capacity());

         system().setUnitCell(parameters);
      }
   }

   /*
   * Output relevant system details to the iteration log file.
   */
   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && verbose() > 1) {
         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nMesh = system().domain().mesh().size();

         // Transfer stress residuals from device to host
         HostDArray<cudaReal> tempH(nFlexibleParams());
         tempH.copySlice(residual(), nMonomer*nMesh);

         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               double str = tempH[counter] / (-1.0 * scaleStress_);
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

   // Private member functions specific to this implementation

   // Calculate the average value of an array.
   template<int D>
   cudaReal AmIteratorGrid<D>::findAverage(VectorT const & field)
   {  return Reduce::sum(field) / (double) field.capacity(); }

}
}
#endif
