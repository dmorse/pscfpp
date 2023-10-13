#ifndef PSPG_AM_ITERATOR_GRID_TPP
#define PSPG_AM_ITERATOR_GRID_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmIteratorGrid.h"
#include <pspg/System.h>
#include <prdc/cuda/RField.h>

#include <prdc/crystal/UnitCell.h>

#include <pscf/inter/Interaction.h>

#include <util/global.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Prdc;

   // Constructor
   template <int D>
   AmIteratorGrid<D>::AmIteratorGrid(System<D>& system)
   : Iterator<D>(system)
   {}

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
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cField(), b.cField(), a.capacity());
   }

   template <int D>
   double AmIteratorGrid<D>::dotProduct(FieldCUDA const & a, 
                                        FieldCUDA const& b) 
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cField(), b.cField(), n);
      return product;
   }

   template <int D>
   double AmIteratorGrid<D>::maxAbs(FieldCUDA const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cField(), n);
      return (double)max;
   }

   template <int D>
   void AmIteratorGrid<D>::updateBasis(RingBuffer<FieldCUDA> & basis, 
                                   RingBuffer<FieldCUDA> const & hists)
   {
      // Make sure at least two histories are stored
      UTIL_CHECK(hists.size() >= 2);

      const int n = hists[0].capacity();
      FieldCUDA newbasis;
      newbasis.allocate(n);

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(n, nBlocks, nThreads);

      pointWiseBinarySubtract<<<nBlocks,nThreads>>>
            (hists[0].cField(),hists[1].cField(),newbasis.cField(),n);

      basis.append(newbasis);
   }

   template <int D>
   void AmIteratorGrid<D>::addHistories(FieldCUDA& trial, 
                                    RingBuffer<FieldCUDA> const & basis, 
                                    DArray<double> coeffs, int nHist)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(trial.capacity(), nBlocks, nThreads);

      for (int i = 0; i < nHist; i++) {
         pointWiseAddScale<<<nBlocks, nThreads>>>
               (trial.cField(), basis[i].cField(), -1*coeffs[i], trial.capacity());
      }
   }

   template <int D>
   void AmIteratorGrid<D>::addPredictedError(FieldCUDA& fieldTrial, 
                                         FieldCUDA const & resTrial, double lambda)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cField(), resTrial.cField(), lambda, 
          fieldTrial.capacity());
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

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // Pointer to fields on system
      DArray<RField<D>> const * currSys = &system().w().rgrid();

      // Loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<nBlocks,nThreads>>>(curr.cField() + i*nMesh, 
                                          (*currSys)[i].cField(), nMesh);
      }

      // If flexible unit cell, also store unit cell parameters
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         const FSArray<double,6> currParam = system().unitCell().parameters();
         // convert into a cudaReal array
         cudaReal* temp = new cudaReal[nParam];
         for (int k = 0; k < nParam; k++) 
               temp[k] = (cudaReal)scaleStress_*currParam[k];
         
         // Copy parameters to the end of the curr array
         cudaMemcpy(curr.cField() + nMonomer*nMesh, temp, 
                    nParam*sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp;
      }
   }

   template <int D>
   void AmIteratorGrid<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();
      // Compute stress if done
      if (isFlexible_) {
         system().mixture().computeStress(system().domain().waveList());
      }
   }

   template <int D>
   void AmIteratorGrid<D>::getResidual(FieldCUDA& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);
      
      // Initialize residuals to zero. Kernel will take care of potential
      // additional elements (n vs nMesh).
      assignUniformReal<<<nBlocks, nThreads>>>(resid.cField(), 0, n);

      // Compute SCF residuals
      for (int i = 0; i < nMonomer; i++) {
         int startIdx = i*nMesh;
         for (int j = 0; j < nMonomer; j++) {
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (resid.cField() + startIdx,
                 system().c().rgrid(j).cField(),
                 interaction_.chi(i, j),
                 nMesh);
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (resid.cField() + startIdx,
                 system().w().rgrid(j).cField(),
                 -interaction_.p(i, j),
                 nMesh);
         }
      }

      // If ensemble is not canonical, account for incompressibility. 
      if (!system().mixture().isCanonical()) {
         cudaReal factor = 1/(cudaReal)interaction_.sumChiInverse();
         for (int i = 0; i < nMonomer; ++i) {
            subtractUniform<<<nBlocks, nThreads>>>(resid.cField() + i*nMesh,
                                                   factor, nMesh);
         }
      } else {
         for (int i = 0; i < nMonomer; i++) {
            // Find current average 
            cudaReal average = findAverage(resid.cField()+i*nMesh, nMesh);
            // subtract out average to set residual average to zero
            subtractUniform<<<nBlocks, nThreads>>>(resid.cField() + i*nMesh,
                                                               average, nMesh);
         }
      }

      // If variable unit cell, compute stress residuals
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         cudaReal* stress = new cudaReal[nParam];

         for (int i = 0; i < nParam; i++) {
            stress[i] = (cudaReal)(-1*scaleStress_*system().mixture().stress(i));
         }

         cudaMemcpy(resid.cField()+nMonomer*nMesh, stress, 
                    nParam*sizeof(cudaReal), cudaMemcpyHostToDevice);
      }
   }

   template <int D>
   void AmIteratorGrid<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // If canonical, explicitly set homogeneous field components 
      if (system().mixture().isCanonical()) {
         cudaReal average, wAverage, cAverage;
         for (int i = 0; i < nMonomer; i++) {
            // Find current spatial average
            average = findAverage(newGuess.cField() + i*nMesh, nMesh);
            
            // Subtract average from field, setting average to zero
            subtractUniform<<<nBlocks, nThreads>>>(newGuess.cField() + i*nMesh,
                                                   average, nMesh);
            
            // Compute the new average omega value, add it to all elements
            wAverage = 0;
            for (int j = 0; j < nMonomer; j++) {
               // Find average concentration for j monomers
               cAverage = findAverage(system().c().rgrid(j).cField(), nMesh);
               wAverage += interaction_.chi(i,j) * cAverage;
            }
            addUniform<<<nBlocks, nThreads>>>(newGuess.cField() + i*nMesh, 
                                              wAverage, nMesh); 
         }
      }

      system().setWRGrid(newGuess);
      system().symmetrizeWFields();

      // If flexible unit cell, update cell parameters 
      if (isFlexible_) {
         FSArray<double,6> parameters;
         const int nParam = system().unitCell().nParameter();
         cudaReal* temp = new cudaReal[nParam];

         cudaMemcpy(temp, newGuess.cField() + nMonomer*nMesh, 
                    nParam*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress_ * (double)temp[i]);
         }
         system().setUnitCell(parameters);

         delete[] temp;
      }

   }

   template<int D>
   void AmIteratorGrid<D>::outputToLog()
   {
      if (isFlexible_) {
         const int nParam = system().unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Parameter " << i << " = "
                        << Dbl(system().unitCell().parameters()[i])
                        << "\n";
         }
      }
   }

   // --- Private member functions that are specific to this implementation --- 

   template<int D> 
   cudaReal AmIteratorGrid<D>::findAverage(cudaReal const * field, int n) 
   {
      cudaReal average = gpuSum(field, n)/n;
      return average;
   }



}
}
#endif
