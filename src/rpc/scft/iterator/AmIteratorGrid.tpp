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
#include <prdc/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/interaction/Interaction.h>
#include <pscf/iterator/NanException.h>
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
      ParamComposite::setClassName("AmIteratorGrid");
      Iterator<D>::isSymmetric_ = true;
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
         if (Iterator<D>::nFlexibleParams() == 0) {
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

      // Read option mixing parameters (lambda, useLambdaRamp, r)
      AmIterTmplT::readMixingParameters(in);

      // Allocate local modified copy of Interaction class
      interaction_.setNMonomer(system().mixture().nMonomer());
   }

   /*
   * Output timing results to log file.
   */
   template <int D>
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
         nEle += Iterator<D>::nFlexibleParams();
      }
      return nEle;
   }

   /*
   * Check if the system has an initial field guess.
   */
   template <int D>
   bool AmIteratorGrid<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field vector (w fields and lattice parameters).
   */
   template <int D>
   void AmIteratorGrid<D>::getCurrent(VectorT& state)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const int nEle = nElements();
      UTIL_CHECK(state.capacity() == nEle);

      // Copy w-fields into a linear array
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
         HostArrayT<RealT> paramsTmp(nFlex);
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
         HostArrayT<RealT> stressTmp(nFlex);
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
         HostArrayT<RealT> paramTmp(nFlex);
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
   * Output relevant system details to the iteration log.
   */
   template <int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && AmIterTmplT::verbose() > 1) {
         RealT res, str;
         UnitCell<D> const & unitCell = system().domain().unitCell();
         const int nParam = unitCell.nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nMesh = system().domain().mesh().size();
         const int begin = nMonomer*nMesh;
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               res = AmIterTmplT::residual()[begin + counter];
               str = - 1.0 * res / scaleStress_;
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
