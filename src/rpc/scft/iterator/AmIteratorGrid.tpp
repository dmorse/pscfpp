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
         nEle += nFlexibleParams();
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
         UTIL_CHECK(nFlexibleParams() > 0);
         UnitCell<D> const & unitCell = system().domain().unitCell();
         FSArray<double, 6> const & parameters = unitCell.parameters();
         const int nParam = unitCell.nParameter();
         DArray<double> paramsTmp(nFlexibleParams());
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

      // Local variables

      // Initialize residual vector to zero
      VecOp::eqS(resid, 0.0);

      // Array of VectorT arrays associated with slices of resid
      DArray<VectorT> residSlices;
      residSlices.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
        residSlices[i].associate(resid, i*nMesh, nMesh);
      }

      // Compute SCF residual vector elements
      double chi, p;
      int begin;
      //int k;
      for (int i = 0; i < nMonomer; ++i) {
         //begin = i*nMesh;
         for (int j = 0; j < nMonomer; ++j) {

            RField<D> const & c = system().c().rgrid(j);
            chi = interaction_.chi(i,j);
            VecOp::addEqVc(residSlices[i], c, chi);

            RField<D> const & w = system().w().rgrid(j);
            p = interaction_.p(i,j);
            VecOp::addEqVc(residSlices[i], w, -1.0*p);

            if (system().h().hasData()) {
               RField<D> const & h = system().h().rgrid(j);
               VecOp::addEqVc(residSlices[i], h, p);
            }

            #if 0
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
            #endif
         }
      }

      // Add term proportional to sum of all concentrations
      double shift = -1.0 / interaction_.sumChiInverse();
      if (system().mask().hasData()) {
         RField<D> const & mask = system().mask().rgrid();
         for (int i = 0; i < nMonomer; ++i) {
            VecOp::addEqVc(residSlices[i], mask, shift);
            #if 0
            begin = i*nMesh;
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] += shift*mask[k];
            }
            #endif
         }
      } else {
         for (int i = 0; i < nMonomer; ++i) {
            VecOp::addEqS(residSlices[i], shift);
            #if 0
            begin = i*nMesh;
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] += shift;
            }
            #endif
         }
      }

      // If system is canonical, set all spatial averages to zero
      if (system().mixture().isCanonical()) {
         double average;
         for (int i = 0; i < nMonomer; ++i) {

            average = computeAverage(residSlices[i]);
            VecOp::subEqS(residSlices[i], average);

            #if 0
            begin = i*nMesh;
            average = 0.0;
            for (k = 0; k < nMesh; ++k) {
               average += resid[begin + k];
            }
            average /= double(nMesh);
            for (k = 0; k < nMesh; ++k) {
               resid[begin + k] -= average;
            }
            #endif
         }
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

         const double coeff = -1.0 * scaleStress_;
         const int nParam = system().domain().unitCell().nParameter();
         begin = nMonomer*nMesh;
         int counter = 0;
         for (int i = 0; i < nParam ; i++) {
            if (flexibleParams_[i]) {
               resid[begin + counter] = coeff * Iterator<D>::stress(i);
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
   void AmIteratorGrid<D>::update(VectorT& newState)
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
         VectorT cAve;
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
      if (isFlexible_) {

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
   template <int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && AmIterTmplT::verbose() > 1) {
         double res, str;
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

   // Private member function

   template <int D>
   double AmIteratorGrid<D>::computeAverage(VectorT const & field) const
   {  return Reduce::sum(field) / double(field.capacity()); }

}
}
#endif
