#ifndef PRDC_CUDA_FFT_BATCHED_H
#define PRDC_CUDA_FFT_BATCHED_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include "RFieldDft.h"
#include <pscf/math/IntVec.h>
#include <pscf/cuda/GpuTypes.h>
#include <util/global.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class FFTBatched
   {

   public:

      /**
      * Default constructor.
      */
      FFTBatched();

      /**
      * Destructor.
      */
      virtual ~FFTBatched();

      /**
      * Set up FFT calculation (get grid dimensions and make FFT plan)
      *
      * \param rDim  dimensions of r-space grid
      * \param kDim  dimensions of k-space grid
      * \param batchSize  number of simultaneous FFTs to perform
      */
      void setup(const IntVec<D>& rDim, const IntVec<D>& kDim, 
                 int batchSize);

      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-grid (device mem)
      * \param out  array of complex values on k-grid (device mem)
      */
      void forwardTransform(DeviceArray<cudaReal>& in, 
                            DeviceArray<cudaComplex>& out) const;

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  ptr to array of complex values on k-grid (device mem)
      * \param out  ptr to array of real values on r-grid (device mem)
      */
      void inverseTransform(DeviceArray<cudaComplex>& in, 
                            DeviceArray<cudaReal>& out) const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      * Set the batch size to a new value. isSetup() must already be true.
      * 
      * \param batchSize  The new batch size.
      */
      void setBatchSize(int batchSize);

      /**
      * Has the setup method been called?
      */
      bool isSetup() const;

   private:

      /// Vector containing number of r-grid points in each direction.
      IntVec<D> meshDimensions_;

      /// Vector containing number of k-grid points in each direction.
      IntVec<D> kMeshDimensions_;

      /// Number of FFTs in batch
      int batchSize_;

      /// Number of points in r-space grid
      int rSize_;

      /// Number of points in k-space grid
      int kSize_;

      /// Pointer to a plan for a forward transform.
      cufftHandle fPlan_;

      /// Pointer to a plan for an inverse transform.
      cufftHandle iPlan_;

      /// Have array dimension and plan been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans for transform and inverse transform.
      * 
      * \param rDim  dimensions of r-space grid
      * \param kDim  dimensions of k-space grid
      * \param batchSize  number of simultaneous FFTs to perform
      */
      void makePlans(const IntVec<D>& rDim, const IntVec<D>& kDim, 
                     int batchSize);

   };

   // Return the dimensions of the grid for which this was allocated.
   template <int D>
   inline const IntVec<D>& FFTBatched<D>::meshDimensions() const
   {  return meshDimensions_; }

   // Has the setup method been called?
   template <int D>
   inline bool FFTBatched<D>::isSetup() const
   { return isSetup_; }

   #ifndef PRDC_CUDA_FFT_BATCHED_TPP
   // Suppress implicit instantiation
   extern template class FFTBatched<1>;
   extern template class FFTBatched<2>;
   extern template class FFTBatched<3>;
   #endif

}
}
}

#endif
