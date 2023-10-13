#ifndef PSPG_FFT_BATCHED_H
#define PSPG_FFT_BATCHED_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/gpu/RField.h>
#include <prdc/gpu/RFieldDft.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <pscf/cuda/GpuResources.h>

//temporary for debugging
#include <iostream>

namespace Pscf {
namespace Pspg {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Pspg_Field_Module
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
      * Check and setup grid dimensions if necessary.
      *
      * \param rField real data on r-space grid (device mem)
      * \param kField complex data on k-space grid (device mem)
      */
      void setup(RField<D>& rField, RFieldDft<D>& kField);

      void setup(const IntVec<D>& rDim,const IntVec<D>& kDim, int batchSize);
      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-space grid (device mem)
      * \param out  array of complex values on k-space grid (device mem)
      */
      void forwardTransform(RField<D>& in, RFieldDft<D>& out);

      void forwardTransform(cudaReal* in, cudaComplex* out, int batchSize);
      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid (device mem)
      * \param out  array of real values on r-space grid (device mem)
      */
      void inverseTransform(RFieldDft<D>& in, RField<D>& out);

      void inverseTransform(cudaComplex* in, cudaReal* out, int batchSize);

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      bool isSetup() const;

      cufftHandle& fPlan();

      cufftHandle& iPlan();


   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid
      int kSize_;

      // Pointer to a plan for a forward transform.
      cufftHandle fPlan_;

      // Pointer to a plan for an inverse transform.
      cufftHandle iPlan_;

      // Have array dimension and plan been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans for transform and inverse transform.
      */
      void makePlans(RField<D>& rField, RFieldDft<D>& kField);
      void makePlans(const IntVec<D>& rDim, const IntVec<D>& kField, int batchSize);

   };

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline const IntVec<D>& FFTBatched<D>::meshDimensions() const
   {  return meshDimensions_; }

   template <int D>
   inline bool FFTBatched<D>::isSetup() const
   { return isSetup_; }

   template <int D>
   inline cufftHandle& FFTBatched<D>::fPlan() 
   { return fPlan_; }

   template <int D>
   inline cufftHandle& FFTBatched<D>::iPlan()
   { return iPlan_; }

   #ifndef PSPG_FFT_BATCHED_TPP
   // Suppress implicit instantiation
   extern template class FFTBatched<1>;
   extern template class FFTBatched<2>;
   extern template class FFTBatched<3>;
   #endif
}
}
//#include "FFTBatched.tpp"
#endif
