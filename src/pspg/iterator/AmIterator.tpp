#ifndef PSPG_AM_ITERATOR_TPP
#define PSPG_AM_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "AmIterator.h"
#include <pspg/field/RDField.h>
#include <pspg/System.h>
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

   template <int D>
   AmIterator<D>::AmIterator(System<D>& system)
   : Iterator<D>(system)
   {}

   template <int D>
   AmIterator<D>::~AmIterator()
   {}

   template <int D>
   double AmIterator<D>::findNorm(FieldCUDA const & hist) 
   {
      const int n = hist.capacity();
      double normResSq = (double)gpuInnerProduct(hist.cDField(), hist.cDField(), n);

      return sqrt(normResSq);
   }

   template <int D>
   double AmIterator<D>::findMaxAbs(FieldCUDA const & hist)
   {
      // use parallel reduction to find maximum.

      // number of data points, each step of the way.
      int n = hist.capacity();
      cudaReal max = gpuMaxAbs(hist.cDField(), n);

      return (double)max;

   }

   template <int D>
   void AmIterator<D>::updateBasis(RingBuffer<FieldCUDA> & basis, RingBuffer<FieldCUDA> const & hists)
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
   double AmIterator<D>::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, int n, int m)
   {      
      return (double)gpuInnerProduct(resBasis[n].cDField(),resBasis[m].cDField(), resBasis[n].capacity());
   }

   template <int D>
   double AmIterator<D>::computeVDotProd(FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int m)
   {
      return (double)gpuInnerProduct(resCurrent.cDField(), resBasis[m].cDField(), resCurrent.capacity());
   }

   template <int D>
   void AmIterator<D>::updateU(DMatrix<double> & U, RingBuffer<FieldCUDA> const & resBasis, int nHist)
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
   void AmIterator<D>::updateV(DArray<double> & v, FieldCUDA const & resCurrent, RingBuffer<FieldCUDA> const & resBasis, int nHist)
   {
      // Compute U matrix's new row 0 and col 0
      // Also, compute each element of v_ vector
      for (int m = 0; m < nHist; ++m) {
         v[m] = computeVDotProd(resCurrent,resBasis,m);
      }
   }

   template <int D>
   void AmIterator<D>::setEqual(FieldCUDA& a, FieldCUDA const & b)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(a.capacity(), nBlocks, nThreads);
      
      UTIL_CHECK(b.capacity() == a.capacity());
      assignReal<<<nBlocks, nThreads>>>(a.cDField(), b.cDField(), a.capacity());
   }

   template <int D>
   void AmIterator<D>::addHistories(FieldCUDA& trial, RingBuffer<FieldCUDA> const & basis, DArray<double> coeffs, int nHist)
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
   void AmIterator<D>::addPredictedError(FieldCUDA& fieldTrial, FieldCUDA const & resTrial, double lambda)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cDField(), resTrial.cDField(), lambda, fieldTrial.capacity());
   }

   template <int D>
   bool AmIterator<D>::hasInitialGuess()
   { return sys_->hasWFields(); }
   
   template <int D>
   int AmIterator<D>::nElements()
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      int nEle = nMonomer*nMesh;

      if (sys_->domain().isFlexible()) {
         nEle += sys_->unitCell().nParameter();
      }

      return nEle;
   }

   template <int D>
   void AmIterator<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();
      const int n = nElements();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // pointer to fields on system
      DArray<RDField<D>> * const currSys = &sys_->wFieldsRGrid();

      // loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<nBlocks,nThreads>>>(curr.cDField() + i*nMesh, (*currSys)[i].cDField(), nMesh);
      }

      // if flexible unit cell, also store unit cell parameters
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         const FSArray<double,6> currParam = sys_->unitCell().parameters();
         // convert into a cudaReal array
         cudaReal* temp = new cudaReal[nParam];
         for (int k = 0; k < nParam; k++) 
               temp[k] = (cudaReal)scaleStress*currParam[k];
         
         // copy paramters to the end of the curr array
         cudaMemcpy(curr.cDField() + nMonomer*nMesh, temp, nParam*sizeof(cudaReal), cudaMemcpyHostToDevice);
         delete[] temp;
      }
   }

   template <int D>
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      sys_->compute();
      // Compute stress if done
      if (sys_->domain().isFlexible()) {
         sys_->mixture().computeStress(sys_->wavelist());
      }
   }

   template <int D>
   void AmIterator<D>::getResidual(FieldCUDA& resid)
   {
      const int n = nElements();
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

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
            pointWiseAddScale<<<nBlocks, nThreads>>>(resid.cDField() + startIdx,
                                                      sys_->cFieldRGrid(j).cDField(),
                                                      sys_->interaction().chi(i, j),
                                                      nMesh);
            pointWiseAddScale<<<nBlocks, nThreads>>>(resid.cDField() + startIdx,
                                                      sys_->wFieldRGrid(j).cDField(),
                                                      -sys_->interaction().idemp(i, j),
                                                      nMesh);
         }
      }
      

      // If not canonical, account for incompressibility. 
      if (!sys_->mixture().isCanonical()) {
         cudaReal factor = 1/(cudaReal)sys_->interaction().sum_inv();
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
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         cudaReal* stress = new cudaReal[nParam];

         for (int i = 0; i < nParam; i++) {
            stress[i] = (cudaReal)(-1*scaleStress*sys_->mixture().stress(i));
         }

         cudaMemcpy(resid.cDField()+nMonomer*nMesh, stress, nParam*sizeof(cudaReal), cudaMemcpyHostToDevice);
      }
   }

   template <int D>
   void AmIterator<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // pointer to field objects in system
      DArray<RDField<D>> * sysfields = &sys_->wFieldsRGrid();

      // Manually and explicitly set homogeneous components of field if canonical
      if (sys_->mixture().isCanonical()) {
         cudaReal average, wAverage, cAverage;
         for (int i = 0; i < nMonomer; i++) {
            // Find current spatial average
            average = findAverage(newGuess.cDField() + i*nMesh, nMesh);
            
            // Subtract average from field, setting average to zero
            subtractUniform<<<nBlocks, nThreads>>>(newGuess.cDField() + i*nMesh, average, nMesh);
            
            // Compute the new average omega value, add it to all elements
            wAverage = 0;
            for (int j = 0; j < nMonomer; j++) {
               // Find average concentration
               cAverage = findAverage(sys_->cFieldRGrid(j).cDField(), nMesh);
               wAverage += sys_->interaction().chi(i,j) * cAverage;
            }
            addUniform<<<nBlocks, nThreads>>>(newGuess.cDField() + i*nMesh, wAverage, nMesh); 
         }
      }

      // copy over grid points
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<nBlocks, nThreads>>>
               ((*sysfields)[i].cDField(), newGuess.cDField() + i*nMesh, nMesh);
      }

      // if flexible unit cell, update parameters well
      if (sys_->domain().isFlexible()) {
         FSArray<double,6> parameters;
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         cudaReal* temp = new cudaReal[nParam];

         cudaMemcpy(temp, newGuess.cDField() + nMonomer*nMesh, nParam*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress * (double)temp[i]);
         }

         sys_->unitCell().setParameters(parameters);            
         sys_->mixture().setupUnitCell(sys_->unitCell(), sys_->wavelist());
         sys_->wavelist().computedKSq(sys_->unitCell());

         delete[] temp;
      }

      // TEMPORARY. PASS THE INPUTS THROUGH A BASIS FILTER.
      sys_->fieldIo().convertRGridToBasis(sys_->wFieldsRGrid(), sys_->wFields());
      sys_->fieldIo().convertRGridToBasis(sys_->cFieldsRGrid(), sys_->cFields());
      sys_->fieldIo().convertBasisToRGrid(sys_->wFields(), sys_->wFieldsRGrid());
      sys_->fieldIo().convertBasisToRGrid(sys_->cFields(), sys_->cFieldsRGrid());

   }

   template<int D>
   void AmIterator<D>::outputToLog()
   {
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         for (int i = 0; i < nParam; i++) {
            Log::file() << "Parameter " << i << " = "
                        << Dbl(sys_->unitCell().parameters()[i])
                        << "\n";
         }
      }
   }

   // --- Private member functions that are specific to this implementation --- 

   template<int D> 
   cudaReal AmIterator<D>::findAverage(cudaReal * const field, int n) 
   {
      cudaReal average = gpuSum(field, n)/n;

      return average;
   }



}
}
#endif
