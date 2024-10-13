#ifndef PRDC_CUDA_FFT_H
#define PRDC_CUDA_FFT_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>
#include <prdc/cuda/CField.h>
#include <prdc/cuda/RFieldDft.h>

//#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class FFT 
   {

   public:

      /**
      * Default constructor.
      */
      FFT();

      /**
      * Destructor.
      */
      virtual ~FFT();

      /**
      * Setup grid dimensions, plans and work space.
      *
      * \param meshDimensions Dimensions of real-space grid.
      */
      void setup(IntVec<D> const & meshDimensions);

      /**
      * Compute forward (real-to-complex) discrete Fourier transform.
      *
      * \param rField  real values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransform(RField<D> & rField, RFieldDft<D>& kField) 
      const;

      /**
      * Compute forward Fourier transform without destroying input.
      *
      * \param rField  real values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransformSafe(RField<D> const & rField, 
                                RFieldDft<D>& kField) const;

      /**
      * Compute inverse (complex-to-real) discrete Fourier transform.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransform(RFieldDft<D> & kField, RField<D>& rField) 
      const;

      /**
      * Compute inverse (complex to real) DFT without destroying input.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransformSafe(RFieldDft<D> const & kField, 
                                RField<D>& rField) const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      *  Has this FFT object been setup?  
      */
      bool isSetup() const;

      #if 0
      /**
      *  Get the plan for the forward DFT.
      */
      cufftHandle& fPlan();

      /**
      *  Get the plan for the inverse DFT.
      */
      cufftHandle& iPlan();
      #endif

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Private RField<D> real array for work space
      mutable RField<D> rFieldCopy_;

      // Private CField<D> real array for work space
      mutable CField<D> cFieldCopy_;

      // Private RFieldDft<D> k-space array work space
      mutable RFieldDft<D> kFieldCopy_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid for transform of real data
      int kSize_;

      // Pointer to a plan for a real-to-complex forward transform
      cufftHandle rcfPlan_;

      // Pointer to a plan for a complex-to-real inverse transform
      cufftHandle criPlan_;

      // Pointer to a plan for a complex-to-complex transform
      cufftHandle ccPlan_;

      // Have array dimension and plans been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans all transform types.
      */
      void makePlans();

   };

   // Declarations of explicit specializations

   template <>
   void FFT<1>::makePlans();

   template <>
   void FFT<2>::makePlans();

   template <>
   void FFT<3>::makePlans();

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline IntVec<D> const & FFT<D>::meshDimensions() const
   {  return meshDimensions_; }

   template <int D>
   inline bool FFT<D>::isSetup() const
   { return isSetup_; }

   #if 0
   template <int D>
   inline cufftHandle& FFT<D>::fPlan() 
   { return rcfPlan_; }

   template <int D>
   inline cufftHandle& FFT<D>::iPlan()
   { return criPlan_; }
   #endif

   #ifndef PRDC_CUDA_FFT_TPP
   // Suppress implicit instantiation
   extern template class FFT<1>;
   extern template class FFT<2>;
   extern template class FFT<3>;
   #endif

   static __global__ 
   void scaleRealData(cudaReal* data, cudaReal scale, int size) {
      //write code that will scale
      int nThreads = blockDim.x * gridDim.x;
      int startId = blockIdx.x * blockDim.x + threadIdx.x;
      for(int i = startId; i < size; i += nThreads ) {
         data[i] *= scale;
      }
   }

}
}
}
#endif
