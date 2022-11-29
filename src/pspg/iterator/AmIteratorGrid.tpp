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
#include <pscf/inter/Interaction.h>
#include <pspg/field/RDField.h>
#include <util/global.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

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

      // Default parameter values
      isFlexible_ = 0;
      scaleStress_ = 10.0;

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
   }

   // Virtual functions used to implement AM algorithm

   template <int D>
   void AmIteratorGrid<D>::setEqual(FieldCUDA& a, FieldCUDA const & b)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cDField(), b.cDField(), a.capacity());
   }

   template <int D>
   double AmIteratorGrid<D>::dotProduct(FieldCUDA const & a, 
                                        FieldCUDA const& b) 
   {
      const int n = a.capacity();
      UTIL_CHECK(b.capacity() == n);
      double product = (double)gpuInnerProduct(a.cDField(), b.cDField(), n);
      return product;
   }

   template <int D>
   double AmIteratorGrid<D>::maxAbs(FieldCUDA const & a)
   {
      int n = a.capacity();
      cudaReal max = gpuMaxAbs(a.cDField(), n);
      return (double)max;
   }

   #if 0
   template <int D>
   double AmIteratorGrid<D>::norm(FieldCUDA const & a) 
   {
      const int n = a.capacity();
      double normResSq;
      normResSq = (double)gpuInnerProduct(a.cDField(), a.cDField(), n);
      return sqrt(normResSq);
   }

   template <int D>
   double 
   AmIteratorGrid<D>::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, 
                                  int n, int m)
   {
      const int n = resBasis[n].capacity();
      UTIL_CHECK(resBasis[m].capacity() == n);
      return (double)gpuInnerProduct(resBasis[n].cDField(),
                                     resBasis[m].cDField(), 
                                     resBasis[n].capacity());
   }

   template <int D>
   double AmIteratorGrid<D>::computeVDotProd(FieldCUDA const & resCurrent, 
                                       RingBuffer<FieldCUDA> const & resBasis, 
                                       int m)
   {
      return (double)gpuInnerProduct(resCurrent.cDField(), 
                                     resBasis[m].cDField(), 
                                     resCurrent.capacity());
   }

   template <int D>
   void AmIteratorGrid<D>::updateU(DMatrix<double> & U, 
                               RingBuffer<FieldCUDA> const & resBasis, 
                               int nHist)
   {
      // Update matrix U by shifting elements diagonally
      int maxHist = U.capacity1();
      for (int m = maxHist-1; m > 0; --m) {
         for (int n = maxHist-1; n > 0; --n) {
            U(m,n) = U(m-1,n-1); 
         }
      }

      // Compute U matrix's new row 0 and col 0
      for (int m = 0; m < nHist; ++m) {
         double dotprod = computeUDotProd(resBasis,m,0);
         U(m,0) = dotprod;
         U(0,m) = dotprod;
      }
   }

   template <int D>
   void AmIteratorGrid<D>::updateV(DArray<double> & v, 
                               FieldCUDA const & resCurrent, 
                               RingBuffer<FieldCUDA> const & resBasis, 
                               int nHist)
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }
   #endif

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
            (hists[0].cDField(),hists[1].cDField(),newbasis.cDField(),n);

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
               (trial.cDField(), basis[i].cDField(), -1*coeffs[i], trial.capacity());
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
         (fieldTrial.cDField(), resTrial.cDField(), lambda, 
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
      DArray<RDField<D>> const * currSys = &system().w().rgrid();

      // Loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<nBlocks,nThreads>>>(curr.cDField() + i*nMesh, 
                                          (*currSys)[i].cDField(), nMesh);
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
         cudaMemcpy(curr.cDField() + nMonomer*nMesh, temp, 
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
      assignUniformReal<<<nBlocks, nThreads>>>(resid.cDField(), 0, n);

      // Compute SCF residuals
      for (int i = 0; i < nMonomer; i++) {
         int startIdx = i*nMesh;
         for (int j = 0; j < nMonomer; j++) {
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (resid.cDField() + startIdx,
                 system().c().rgrid(j).cDField(),
                 system().interaction().chi(i, j),
                 nMesh);
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (resid.cDField() + startIdx,
                 system().w().rgrid(j).cDField(),
                 -system().interaction().idemp(i, j),
                 nMesh);
         }
      }

      // If ensemble is not canonical, account for incompressibility. 
      if (!system().mixture().isCanonical()) {
         cudaReal factor = 1/(cudaReal)system().interaction().sum_inv();
         for (int i = 0; i < nMonomer; ++i) {
            subtractUniform<<<nBlocks, nThreads>>>(resid.cDField() + i*nMesh,
                                                   factor, nMesh);
         }
      } else {
         for (int i = 0; i < nMonomer; i++) {
            // Find current average 
            cudaReal average = findAverage(resid.cDField()+i*nMesh, nMesh);
            // subtract out average to set residual average to zero
            subtractUniform<<<nBlocks, nThreads>>>(resid.cDField() + i*nMesh,
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

         cudaMemcpy(resid.cDField()+nMonomer*nMesh, stress, 
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
            average = findAverage(newGuess.cDField() + i*nMesh, nMesh);
            
            // Subtract average from field, setting average to zero
            subtractUniform<<<nBlocks, nThreads>>>(newGuess.cDField() + i*nMesh,
                                                   average, nMesh);
            
            // Compute the new average omega value, add it to all elements
            wAverage = 0;
            for (int j = 0; j < nMonomer; j++) {
               // Find average concentration for j monomers
               cAverage = findAverage(system().c().rgrid(j).cDField(), nMesh);
               wAverage += system().interaction().chi(i,j) * cAverage;
            }
            addUniform<<<nBlocks, nThreads>>>(newGuess.cDField() + i*nMesh, 
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

         cudaMemcpy(temp, newGuess.cDField() + nMonomer*nMesh, 
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
