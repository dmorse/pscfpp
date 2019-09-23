#ifndef PSSP_GPU_FFT_BATCHED_H
#define PSSP_GPU_FFT_BATCHED_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp_gpu/field/RDField.h>
#include <pssp_gpu/field/RDFieldDft.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>

//temporary for debugging
#include <iostream>
#include <pssp_gpu/field/FFT.h> //for definition of rtype
//#ifdef SINGLE_PRECISION
//typedef float rtype;
//#else
//typedef double rtype;
//#endif
//typedef rtype; //some type in FFT.h

namespace Pscf {
namespace Pssp_gpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Pssp_Field_Module
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
      * \param rDField real data on r-space grid (device mem)
      * \param kDField complex data on k-space grid (device mem)
      */
      void setup(RDField<D>& rDField, RDFieldDft<D>& kDField);

      void setup(const IntVec<D>& rDim,const IntVec<D>& kDim, int batchSize);
      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-space grid (device mem)
      * \param out  array of complex values on k-space grid (device mem)
      */
      void forwardTransform(RDField<D>& in, RDFieldDft<D>& out);

      void forwardTransform(cufftReal* in, cufftComplex* out, int batchSize);
      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid (device mem)
      * \param out  array of real values on r-space grid (device mem)
      */
      void inverseTransform(RDFieldDft<D>& in, RDField<D>& out);

      void inverseTransform(cufftComplex* in, cufftReal* out, int batchSize);

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
      void makePlans(RDField<D>& rDField, RDFieldDft<D>& kDField);
      void makePlans(const IntVec<D>& rDim, const IntVec<D>& kDField, int batchSize);

   };

   // Declarations of explicit specializations

   //template <>
   //void FFT<1>::makePlans(RDField<1>& rField, RDFieldDft<1>& kField);

   //template <>
   //void FFT<2>::makePlans(RDField<2>& rField, RDFieldDft<2>& kField);

   template <>
   void FFTBatched<3>::makePlans(RDField<3>& rField, RDFieldDft<3>& kField);

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

}
}
#include "FFTBatched.tpp"
#endif
