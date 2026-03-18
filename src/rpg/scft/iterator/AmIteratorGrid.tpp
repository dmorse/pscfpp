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
   * Get the current state vector (w fields and lattice parameters).
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
         int nFlex = Iterator<D>::nFlexibleParams();
         UTIL_CHECK(nFlex > 0);
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double, 6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         DArray<cudaReal> paramsTmp(nFlex);
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               paramsTmp[counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlex);

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
      // Precondition 
      const int n = nElements();
      UTIL_CHECK(resid.capacity() == n);

      // Constants
      const RealT shift = -1.0 / interaction_.sumChiInverse();
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const bool hasHext = system().h().hasData();
      const bool hasMask = system().mask().hasData();
      const bool isCanonical = system().mixture().isCanonical();

      // Initialize residual vector to zero
      VecOp::eqS(resid, 0.0);

      // Compute field residual elements
      VectorT slice;
      RealT chi, p, average;
      for (int i = 0; i < nMonomer; ++i) {
         slice.associate(resid, i*nMesh, nMesh);

         // Matrix products
         for (int j = 0; j < nMonomer; ++j) {
            chi = interaction_.chi(i,j);
            p = interaction_.p(i,j);
            VecOp::addEqVc(slice, system().c().rgrid(j), chi);
            VecOp::addEqVc(slice, system().w().rgrid(j), -1.0*p);
            if (hasHext) {
               VecOp::addEqVc(slice, system().h().rgrid(j), p);
            }
         }

         // Term proportional to required sum of concentrations
         if (hasMask) {
            VecOp::addEqVc(slice, system().mask().rgrid(), shift);
         } else {
            if (!isCanonical) {
               VecOp::addEqS(slice, shift);
            }
         }

         // If canonical ensemble, subtract average from residual slice
         if (isCanonical) {
            average = Reduce::sum(slice);
            average /= RealT(nMesh);
            VecOp::subEqS(slice, average);
         }

         slice.dissociate();
      }

      // If flexible unit cell, then compute stress residuals
      if (isFlexible_) {

         // Combined -1 factor and stress scaling here. This is okay:
         // - residuals only show up as dot products (U, v, norm)
         //   or with their absolute value taken (max), so the
         //   sign on a given residual vector element is not relevant
         //   as long as it is consistent across all vectors
         // - The scaling is applied here and to the unit cell param
         //   storage, so that updating is done on the same scale,
         //   and then undone right before passing to the unit cell.

         const RealT scale = -1.0 * scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         const int nFlex = Iterator<D>::nFlexibleParams();
         HostDArray<RealT> stressTmp(nFlex);
         int counter = 0;
         for (int i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               stressTmp[counter] = scale * Iterator<D>::stress(i);
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlex);
         UTIL_CHECK(resid.capacity() == (nMonomer * nMesh) + nFlex);

         // Copy stress residuals to the end of the resid array
         slice.associate(resid, nMonomer * nMesh, nFlex);
         slice = stressTmp; // copy from host to device, for GPU code
         slice.dissociate();
      }

   }

   /*
   * Update the current system field vector.
   */
   template <int D>
   void AmIteratorGrid<D>::update(VectorT& newState)
   {
      // Constants and references to system components
      Domain<D> const & domain = system().domain();
      Mesh<D> const & mesh = domain.mesh();
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = mesh.size();

      // Allocate wFields container
      DArray< RField<D> > wFields;
      wFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; i++) {
         wFields[i].allocate(mesh.dimensions());
         UTIL_CHECK(wFields[i].capacity() == nMesh);
      }

      // Copy new fields from newState vector to wFields container
      VectorT slice;
      for (int i = 0; i < nMonomer; i++) {
         slice.associate(newState, i*nMesh, nMesh);
         VecOp::eqV(wFields[i], slice);
         slice.dissociate();
      }

      // If canonical, explicitly set homogeneous field components
      if (system().mixture().isCanonical()) {

         // Subtract spatial average from each w field
         RealT wAve;
         for (int i = 0; i < nMonomer; ++i) {
            wAve = Reduce::sum(wFields[i]);
            wAve /= RealT(nMesh);
            VecOp::subEqS(wFields[i], wAve);
         }

         // Compute spatial averages of all concentration fields
         DArray<RealT> cAve;
         cAve.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            cAve[i] = Reduce::sum(system().c().rgrid(i));
            cAve[i] /= RealT(nMesh);
         }

         // Add average values arising from interactions
         RealT chi;
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
            RealT hAve;
            for (int i = 0; i < nMonomer; ++i) {
               hAve = Reduce::sum(system().h().rgrid(i));
               hAve /= RealT(nMesh);
               VecOp::addEqS(wFields[i], hAve);
            }
         }
      }

      // Set fields in system w container
      system().w().setRGrid(wFields);

      // If flexible, update unit cell parameters
      if (isFlexible_) {
         const int nParam = domain.unitCell().nParameter();
         const int nFlex = Iterator<D>::nFlexibleParams();

         // Initialize parameters array with current values
         FSArray<double, 6> parameters;
         parameters = domain.unitCell().parameters();

         // Copy parameter entries from newState to a local array
         HostDArray<RealT> paramTmp(nFlex);
         slice.associate(newState, nMonomer*nMesh, nFlex);
         paramTmp = slice;
         slice.dissociate();

         // Reset any parameters that are flexible
         RealT scale = 1.0 / scaleStress_;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               parameters[i] = scale * paramTmp[counter];
               counter++;
            }
         }
         UTIL_CHECK(counter == nFlex);

         // Set system unit cell parameters
         system().setUnitCell(parameters);
      }

   }

   /*
   * Output relevant system details to the iteration log file.
   */
   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && AmIterTmplT::verbose() > 1) {
         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nFlex = Iterator<D>::nFlexibleParams();
         const int nMonomer = system().mixture().nMonomer();
         const int nMesh = system().domain().mesh().size();

         // Transfer stress residuals from device to host
         HostDArray<cudaReal> stressTmp(nFlex);
         stressTmp.copySlice(residual(), nMonomer*nMesh);

         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               double str = stressTmp[counter] / (-1.0 * scaleStress_);
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

   #if 0
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
            cudaReal average = findAverage(residSlices[i]);
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
   #endif

   #if 0
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
   #endif

}
}
#endif
