#ifndef PRDC_CUDA_FFT_BATCHED_H
#define PRDC_CUDA_FFT_BATCHED_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include "RFieldDft.h"
#include "types.h"
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Batched Fourier transform wrapper for real data.
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
      * \param meshDimensions  dimensions of real-space grid
      * \param batchSize  number of simultaneous FFTs to perform
      */
      void setup(IntVec<D> const & meshDimensions, int batchSize);
      
      /**
      * Compute batched forward (real-to-complex) DFTs.
      * 
      * The resulting complex array is the output of cuFFT's batched 
      * forward transform method scaled by a factor of 1/N (where N is 
      * the number of grid points). The scaling ensures that a round-trip
      * Fourier transform (forward + inverse) reproduces the original
      * input array.
      * 
      * The number of Fourier transforms performed per batch is pre-set
      * by the setup() method.
      * 
      * This transform is inherently "safe", meaning that the input array 
      * will not be modified during the DFT.
      *
      * \param rFields  real values on r-grid (input, gpu mem)
      * \param kFields  complex values on k-grid (output, gpu mem)
      */
      void forwardTransform(DeviceArray<cudaReal> const & rFields, 
                            DeviceArray<cudaComplex>& kFields) const;

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      * 
      * The resulting real array is the unaltered output of cuFFT's inverse 
      * transform function.
      * 
      * This transform is inherently "unsafe", meaning that the input array
      * will be overwritten during the DFT. Accordingly, the input array
      * kField is not const like it is in all the other transform methods.
      *
      * \param kFields  complex values on k-grid (input, gpu mem)
      * \param rFields  real values on r-grid (output, gpu mem)
      */
      void inverseTransformUnsafe(DeviceArray<cudaComplex>& kFields, 
                                  DeviceArray<cudaReal>& rFields) const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      * Set the batch size to a new value. isSetup() must already be true.
      * 
      * \param batchSize  The new batch size.
      */
      void resetBatchSize(int batchSize);

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
      * Make cuFFT plans for transform and inverse transform.
      * 
      * \param batchSize  number of simultaneous FFTs to perform
      */
      void makePlans(int batchSize);

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
