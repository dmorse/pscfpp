#ifndef RPG_AM_ITERATOR_GRID_TPP
#define RPG_AM_ITERATOR_GRID_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpg/System.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/cuda/resources.h>
#include <pscf/inter/Interaction.h>
#include <pscf/iterator/NanException.h>

#include <util/global.h>
#include <cmath>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   // Constructor
   template <int D>
   AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
    : Iterator<D>(system),
      imposedFields_(system)
   {
      setClassName("AmIteratorGrid");
      isSymmetric_ = false; 
   }

   // Destructor
   template <int D>
   AmIteratorGrid<D>::~AmIteratorGrid()
   {}

   // Read parameter file block 
   template <int D>
   void AmIteratorGrid<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters 
      AmIteratorTmpl<Iterator<D>,FieldCUDA>::readParameters(in);
      AmIteratorTmpl<Iterator<D>,FieldCUDA>::readErrorType(in);

      // Allocate local modified copy of Interaction class
      interaction_.setNMonomer(system().mixture().nMonomer());

      // Default parameter values
      isFlexible_ = 1;
      scaleStress_ = 10.0;

      int np = system().domain().unitCell().nParameter();
      UTIL_CHECK(np > 0);
      UTIL_CHECK(np <= 6);
      UTIL_CHECK(system().domain().unitCell().lattice() != UnitCell<D>::Null);

      // Read in optional isFlexible value
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
         if (nFlexibleParams() == 0) isFlexible_ = false;
      } else { // isFlexible_ = false
         flexibleParams_.clear();
         for (int i = 0; i < np; i++) {
            flexibleParams_.append(false); // Set all values to false
         }
      }

      // Read optional scaleStress value
      readOptional(in, "scaleStress", scaleStress_);

      // Read optional ImposedFieldsGenerator object
      readParamCompositeOptional(in, imposedFields_);
   }

   // Output timing results to log file.
   template<int D>
   void AmIteratorGrid<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Iterator times contributions:\n";
      AmIteratorTmpl<Iterator<D>, FieldCUDA >::outputTimers(out);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorGrid<D>::setup(bool isContinuation)
   {
      if (imposedFields_.isActive()) {
         imposedFields_.setup();
      }

      AmIteratorTmpl<Iterator<D>, FieldCUDA>::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions used to implement AM algorithm

   // Set vector a equal to vector b (a = b).
   template <int D>
   void AmIteratorGrid<D>::setEqual(FieldCUDA& a, FieldCUDA const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   // Compute and return inner product of two real fields.
   template <int D>
   double AmIteratorGrid<D>::dotProduct(FieldCUDA const & a, 
                                        FieldCUDA const & b) 
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      double val = Reduce::innerProduct(a, b);
      if (std::isnan(val)) { // if value is NaN, throw NanException
         throw NanException("AmIteratorGrid::dotProduct", __FILE__, 
                            __LINE__, 0);
      }
      return val;
   }

   // Find the maximum magnitude element of a residual vector.
   template <int D>
   double AmIteratorGrid<D>::maxAbs(FieldCUDA const & a)
   {
      double val = Reduce::maxAbs(a);
      if (std::isnan(val)) { // if value is NaN, throw NanException
         throw NanException("AmIteratorGrid::maxAbs", __FILE__, __LINE__, 0);
      }
      return val;
   }

   // Update the series of residual vectors.
   template <int D>
   void AmIteratorGrid<D>::updateBasis(RingBuffer<FieldCUDA> & basis, 
                                       RingBuffer<FieldCUDA> const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      // Set up array to store new basis
      basis.advance();
      if (basis[0].isAllocated()) {
         UTIL_CHECK(basis[0].capacity() == hists[0].capacity());
      } else {
         basis[0].allocate(hists[0].capacity());
      }

      // New basis vector is difference between two most recent states
      VecOp::subVV(basis[0], hists[0], hists[1]);
   }

   // Compute trial field so as to minimize L2 norm of residual.
   template <int D>
   void 
   AmIteratorGrid<D>::addHistories(FieldCUDA& trial, 
                                   RingBuffer<FieldCUDA> const & basis, 
                                   DArray<double> coeffs, int nHist)
   {
      for (int i = 0; i < nHist; i++) {
         VecOp::addEqVc(trial, basis[i], -1.0 * coeffs[i]);
      }
   }

   // Add predicted error to the trial field.
   template <int D>
   void AmIteratorGrid<D>::addPredictedError(FieldCUDA& fieldTrial, 
                                             FieldCUDA const & resTrial, 
                                             double lambda)
   {
      VecOp::addEqVc(fieldTrial, resTrial, lambda);
   }

   // Checks if the system has an initial guess
   template <int D>
   bool AmIteratorGrid<D>::hasInitialGuess()
   { return system().w().hasData(); }

   // Compute the number of elements in the residual vector.
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

   // Get the current w fields and lattice parameters.
   template <int D>
   void AmIteratorGrid<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();
      const int n = nElements();
      UTIL_CHECK(curr.capacity() == n);

      // Loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         VecOp::eqV(curr, system().w().rgrid(i), i*nMesh, 0, nMesh);
      }

      // If flexible unit cell, also store unit cell parameters
      if (isFlexible_) {
         const int nParam = system().domain().unitCell().nParameter();
         FSArray<double,6> const & parameters
                              = system().domain().unitCell().parameters();

         // convert into a cudaReal array
         HostDArray<cudaReal> tempH(nFlexibleParams());
         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               tempH[counter] = scaleStress_ * parameters[i];
               counter++;
            }
         }
         UTIL_CHECK(counter == tempH.capacity());
         
         // Copy parameters to the end of the curr array
         FieldCUDA tempD;
         tempD.associate(curr, nMonomer*nMesh, tempH.capacity());
         tempD = tempH; // copy from host to device
      }
   }

   // Solve MDE for current state of system.
   template <int D>
   void AmIteratorGrid<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute(isFlexible_);
   }

   // Gets the residual vector from system.
   template <int D>
   void AmIteratorGrid<D>::getResidual(FieldCUDA& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();

      // Initialize residuals to zero. Kernel will take care of potential
      // additional elements (n vs nMesh).
      VecOp::eqS(resid, 0.0);

      // Array of FieldCUDA arrays associated with slices of resid.
      // one FieldCUDA array per monomer species, each of size nMesh.
      DArray<FieldCUDA> residSlices;
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
      if (system().hasMask()) {
         double coeff = -1.0 / interaction_.sumChiInverse();
         for (int i = 0; i < nMonomer; ++i) {
            VecOp::addEqVc(residSlices[i], system().mask().rgrid(), coeff);
         }
      }

      // If iterator has external fields, account for them in the values 
      // of the residuals
      if (system().hasExternalFields()) {
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
         const int nParam = system().domain().unitCell().nParameter();
         HostDArray<cudaReal> stressH(nFlexibleParams());

         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               double stress = system().mixture().stress(i);

               // Correct stress to account for effect of imposed fields
               if (imposedFields_.isActive()) {
                  stress = imposedFields_.correctedStress(i,stress);
               } 

               stressH[counter] = -1 * scaleStress_ * stress;
               counter++;
            }
         }
         UTIL_CHECK(counter == stressH.capacity());

         FieldCUDA stressD;
         stressD.associate(resid, nMonomer*nMesh, stressH.capacity());
         stressD = stressH; // copy from host to device
      }
   }

   // Update the system with a new trial field vector.
   template <int D>
   void AmIteratorGrid<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().domain().mesh().size();

      // If canonical, explicitly set homogeneous field components 
      if (system().mixture().isCanonical()) {
         cudaReal average, wAverage, cAverage;
         for (int i = 0; i < nMonomer; i++) {
            // Define array associated with a slice of newGuess
            FieldCUDA ngSlice;
            ngSlice.associate(newGuess, i*nMesh, nMesh);

            // Find current spatial average
            average = findAverage(ngSlice);
            
            // Subtract average from field, setting average to zero
            VecOp::subEqS(ngSlice, average);
            
            // Compute the new average omega value
            wAverage = 0;
            for (int j = 0; j < nMonomer; j++) {
               // Find average concentration for j monomers
               if (system().w().isSymmetric()) { // c().basis() has data
                  cAverage = system().c().basis(j)[0];
               } else { // average must be calculated
                  cAverage = findAverage(system().c().rgrid(j));
               }
               wAverage += interaction_.chi(i,j) * cAverage;
            }

            // If system has external fields, include them in homogeneous field
            if (system().hasExternalFields()) {
               if (system().h().isSymmetric()) { // h().basis() has data
                  wAverage += system().h().basis(i)[0];
               } else { // average must be calculated
                  wAverage += findAverage(system().h().rgrid(i));
               }
            }

            // Add new average omega value to the field
            VecOp::addEqS(ngSlice, wAverage);
         }
      }

      system().setWRGrid(newGuess);
      system().symmetrizeWFields();

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

      // Update imposed fields if needed
      if (imposedFields_.isActive()) {
         imposedFields_.update();
      }
   }

   // Output relevant system details to the iteration log file.
   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && verbose() > 1) {
         const int nParam = system().domain().unitCell().nParameter();
         const int nMonomer = system().mixture().nMonomer();
         const int nMesh = system().domain().mesh().size();

         // transfer stress residuals from device to host
         HostDArray<cudaReal> tempH(nFlexibleParams());
         tempH.copySlice(residual(), nMonomer*nMesh);

         int counter = 0;
         for (int i = 0; i < nParam; i++) {
            if (flexibleParams_[i]) {
               double stress = tempH[counter] / (-1.0 * scaleStress_);
               Log::file() 
                   << " Cell Param  " << i << " = "
                   << Dbl(system().domain().unitCell().parameters()[i], 15)
                   << " , stress = " 
                   << Dbl(stress, 15)
                   << "\n";
               counter++;
            }
         }
      }
   }

   // Return specialized sweep parameter types to add to the Sweep object
   template<int D>
   GArray<ParameterType> AmIteratorGrid<D>::getParameterTypes()
   {
      GArray<ParameterType> arr;
      if (imposedFields_.isActive()) {
         arr = imposedFields_.getParameterTypes();
      } 
      return arr;
   }
   // Set the value of a specialized sweep parameter
   template<int D>
   void AmIteratorGrid<D>::setParameter(std::string name, DArray<int> ids, 
                                         double value, bool& success)
   {
      if (imposedFields_.isActive()) {
         imposedFields_.setParameter(name, ids, value, success);
      } else {
         success = false;
      }
   }
   // Get the value of a specialized sweep parameter
   template<int D>
   double AmIteratorGrid<D>::getParameter(std::string name, 
                                           DArray<int> ids, bool& success)
   const
   {
      if (imposedFields_.isActive()) {
         return imposedFields_.getParameter(name, ids, success);
      } else {
         success = false;
         return 0.0;
      }
   }

   // --- Private member functions specific to this implementation --- 

   // Calculate the average value of an array.
   template<int D> 
   cudaReal AmIteratorGrid<D>::findAverage(FieldCUDA const & field) 
   {
      return Reduce::sum(field) / field.capacity();
   }

}
}
#endif
