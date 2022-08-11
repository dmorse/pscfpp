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
   void AmIterator<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters 
      AmIteratorTmpl<Iterator<D>,FieldCUDA>::readParameters(in);

      // Default parameter values
      isFlexible_ = 0;
      scaleStress_ = 10.0;

      // Read in additional parameters
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
   }


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
   void AmIterator<D>::updateBasis(RingBuffer<FieldCUDA> & basis, 
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
   double 
   AmIterator<D>::computeUDotProd(RingBuffer<FieldCUDA> const & resBasis, 
                                  int n, int m)
   {
      return (double)gpuInnerProduct(resBasis[n].cDField(),
                                     resBasis[m].cDField(), 
                                     resBasis[n].capacity());
   }

   template <int D>
   double AmIterator<D>::computeVDotProd(FieldCUDA const & resCurrent, 
                                         RingBuffer<FieldCUDA> const & resBasis, 
                                         int m)
   {
      return (double)gpuInnerProduct(resCurrent.cDField(), 
                                     resBasis[m].cDField(), 
                                     resCurrent.capacity());
   }

   template <int D>
   void AmIterator<D>::updateU(DMatrix<double> & U, 
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
   void AmIterator<D>::updateV(DArray<double> & v, 
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
   void AmIterator<D>::addHistories(FieldCUDA& trial, 
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
   void AmIterator<D>::addPredictedError(FieldCUDA& fieldTrial, 
                                         FieldCUDA const & resTrial, double lambda)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(fieldTrial.capacity(), nBlocks, nThreads);

      pointWiseAddScale<<<nBlocks, nThreads>>>
         (fieldTrial.cDField(), resTrial.cDField(), lambda, fieldTrial.capacity());
   }

   template <int D>
   bool AmIterator<D>::hasInitialGuess()
   { return system().hasWFields(); }
   
   template <int D>
   int AmIterator<D>::nElements()
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
   void AmIterator<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();
      const int n = nElements();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // Pointer to fields on system
      DArray<RDField<D>> * const currSys = &system().wFieldsRGrid();

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
   void AmIterator<D>::evaluate()
   {
      // Solve MDEs for current omega field
      system().compute();
      // Compute stress if done
      if (isFlexible_) {
         system().mixture().computeStress(system().wavelist());
      }
   }

   template <int D>
   void AmIterator<D>::getResidual(FieldCUDA& resid)
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
                 system().cFieldRGrid(j).cDField(),
                 system().interaction().chi(i, j),
                 nMesh);
            pointWiseAddScale<<<nBlocks, nThreads>>>
                (resid.cDField() + startIdx,
                 system().wFieldRGrid(j).cDField(),
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
   void AmIterator<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nMesh = system().mesh().size();

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // pointer to field objects in system
      DArray<RDField<D>> * sysfields = &system().wFieldsRGrid();

      // Manually and explicitly set homogeneous components of field if canonical
      if (system().mixture().isCanonical()) {
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
               cAverage = findAverage(system().cFieldRGrid(j).cDField(), nMesh);
               wAverage += system().interaction().chi(i,j) * cAverage;
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
      if (isFlexible_) {
         FSArray<double,6> parameters;
         const int nParam = system().unitCell().nParameter();
         cudaReal* temp = new cudaReal[nParam];

         cudaMemcpy(temp, newGuess.cDField() + nMonomer*nMesh, nParam*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress_ * (double)temp[i]);
         }

         #if 0
         system().unitCell().setParameters(parameters);            
         system().mixture().setupUnitCell(system().unitCell(), system().wavelist());
         system().wavelist().computedKSq(system().unitCell());
         #endif
   
         system().setUnitCell(parameters);

         delete[] temp;
      }

      // TEMPORARY. PASS THE INPUTS THROUGH A BASIS FILTER.
      system().fieldIo().convertRGridToBasis(system().wFieldsRGrid(), system().wFields());
      system().fieldIo().convertRGridToBasis(system().cFieldsRGrid(), system().cFields());
      system().fieldIo().convertBasisToRGrid(system().wFields(), system().wFieldsRGrid());
      system().fieldIo().convertBasisToRGrid(system().cFields(), system().cFieldsRGrid());

   }

   template<int D>
   void AmIterator<D>::outputToLog()
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
   cudaReal AmIterator<D>::findAverage(cudaReal * const field, int n) 
   {
      cudaReal average = gpuSum(field, n)/n;
      return average;
   }



}
}
#endif
