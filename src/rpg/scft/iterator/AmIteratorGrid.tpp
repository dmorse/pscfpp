#ifndef RPG_AM_ITERATOR_GRID_TPP
#define RPG_AM_ITERATOR_GRID_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <rpg/System.h>
#include <prdc/cuda/RField.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/inter/Interaction.h>
#include <pscf/cuda/GpuResources.h>

#include <util/global.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   // Constructor
   template <int D>
   AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
    : Iterator<D>(system)
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

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
   }

   // Protected virtual function

   // Setup before entering iteration loop
   template <int D>
   void AmIteratorGrid<D>::setup(bool isContinuation)
   {
      AmIteratorTmpl<Iterator<D>, FieldCUDA>::setup(isContinuation);
      interaction_.update(system().interaction());
   }

   // Private virtual functions used to implement AM algorithm

   template <int D>
   void AmIteratorGrid<D>::setEqual(FieldCUDA& a, FieldCUDA const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   template <int D>
   double AmIteratorGrid<D>::dotProduct(FieldCUDA const & a, 
                                        FieldCUDA const & b) 
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   template <int D>
   double AmIteratorGrid<D>::maxAbs(FieldCUDA const & a)
   {
      return Reduce::maxAbs(a);
   }

   template <int D>
   void AmIteratorGrid<D>::updateBasis(RingBuffer<FieldCUDA> & basis, 
                                       RingBuffer<FieldCUDA> const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      basis.advance();
      if (basis[0].isAllocated()) {
         UTIL_CHECK(basis[0].capacity() == hists[0].capacity());
      } else {
         basis[0].allocate(hists[0].capacity());
      }

      VecOp::subVV(basis[0], hists[0], hists[1]);
   }

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

   template <int D>
   void AmIteratorGrid<D>::addPredictedError(FieldCUDA& fieldTrial, 
                                             FieldCUDA const & resTrial, 
                                             double lambda)
   {
      VecOp::addEqVc(fieldTrial, resTrial, lambda);
   }

   template <int D>
   bool AmIteratorGrid<D>::hasInitialGuess()
   { return system().w().hasData(); }
   
   template <int D>
   int AmIteratorGrid<D>::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

      int nEle = nMonomer*nMesh;

      if (isFlexible_) {
         nEle += system().unitCell().nParameter();
      }

      return nEle;
   }

   template <int D>
   void AmIteratorGrid<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();
      const int n = nElements();
      UTIL_CHECK(curr.capacity() == n);

      // Loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         VecOp::eqV(curr, system().w().rgrid(i), i*nMesh, 0, nMesh);
      }

      // If flexible unit cell, also store unit cell parameters
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         const FSArray<double, 6> currParam = 
                                       system().unitCell().parameters();

         // convert into a cudaReal array
         HostDArray<cudaReal> tempH(nParam);
         for (int k = 0; k < nParam; k++) {
            tempH[k] = scaleStress_ * currParam[k];
         }
         
         // Copy parameters to the end of the curr array
         FieldCUDA tempD;
         tempD.associate(curr, nMonomer*nMesh, nParam);
         tempD = tempH; // copy from host to device
      }
   }

   template <int D>
   void AmIteratorGrid<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute(isFlexible_);
   }

   template <int D>
   void AmIteratorGrid<D>::getResidual(FieldCUDA& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

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
         const int nParam = system().unitCell().nParameter();
         HostDArray<cudaReal> stressH(nParam);

         for (int i = 0; i < nParam; i++) {
            stressH[i] = (cudaReal)(-1 * scaleStress_ * 
                                   system().mixture().stress(i));
         }

         FieldCUDA stressD;
         stressD.associate(resid, nMonomer*nMesh, nParam);
         stressD = stressH; // copy from host to device
      }
   }

   template <int D>
   void AmIteratorGrid<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

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
            
            // Compute the new average omega value, add it to all elements
            wAverage = 0;
            for (int j = 0; j < nMonomer; j++) {
               // Find average concentration for j monomers
               cAverage = findAverage(system().c().rgrid(j));
               wAverage += interaction_.chi(i,j) * cAverage;
            }
            VecOp::addEqS(ngSlice, wAverage);
         }
      }

      system().setWRGrid(newGuess);
      system().symmetrizeWFields();

      // If flexible unit cell, update cell parameters 
      if (isFlexible_) {
         FSArray<double,6> parameters;
         const int nParam = system().unitCell().nParameter();
         HostDArray<cudaReal> tempH(nParam);
         FieldCUDA tempD;
         tempD.associate(newGuess, nMonomer*nMesh, nParam);
         tempH = tempD; // transfer from device to host
         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress_ * (double)tempH[i]);
         }
         system().setUnitCell(parameters);
      }

   }

   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_ && verbose() > 1) {
         const int nParam = system().unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Parameter " << i << " = "
                        << Dbl(system().unitCell().parameters()[i])
                        << "\n";
         }
      }
   }

   // --- Private member functions specific to this implementation --- 

   template<int D> 
   cudaReal AmIteratorGrid<D>::findAverage(FieldCUDA const & field) 
   {
      return Reduce::sum(field) / field.capacity();
   }

}
}
#endif
