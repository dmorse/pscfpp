#ifndef PSSP_FFT_H
#define PSSP_FFT_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp/field/RField.h>
#include <pssp/field/RFieldDFT.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pssp
{

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform plan for real data.
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
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-space grid
      * \param out array of complex values on k-space grid
      */
      void forwardTransform(RField<D>& in, RFieldDFT<D>& out);

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid
      * \param out  array of real values on r-space grid
      * \param meshDimensions number of grid points in each direction.
      */
      void inverseTransform(RFieldDFT<D>& in, RField<D>& out);

      /**
      * Get the dimensions of the grid for which this was allocated.
      *
      * \throw Exception if dimensions of space do not match.
      *
      * \param dimensions number of grid points in each direction.
      */
      void getMeshDimensions(IntVec<D>& meshDimensions) const;

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid
      int kSize_;

      // Pointer to a plan for a forward transform.
      fftw_plan fPlan_;

      // Pointer to a plan for an inverse transform.
      fftw_plan rPlan_;

      // True when grid dimensions have been set.
      bool hasDimensions_;

      /**
      * Check and setup grid dimensions if necessary.
      *
      * \param dimensions number of grid points in each direction.
      */
      void setup(const RField<D>& rField, const RFieldDFT<D>& kField);

      /**
      * Set new grid dimensions.
      *
      * \param dimensions number of grid points in each direction.
      */
      void setDimensions(const IntVec<D>& meshDimensions);

   };

   // Declarations of explicit specializations

   template <>
   void FFT<1>::forwardTransform(RField<1>& in, RFieldDFT<1>& out);

   template <>
   void FFT<2>::forwardTransform(RField<2>& in, RFieldDFT<2>& out);

   template <>
   void FFT<3>::forwardTransform(RField<3>& in, RFieldDFT<3>& out);

   template <>
   void FFT<1>::inverseTransform(RFieldDFT<1>& in, RField<1>& out);

   template <>
   void FFT<2>::inverseTransform(RFieldDFT<2>& in, RField<2>& out);

   template <>
   void FFT<3>::inverseTransform(RFieldDFT<3>& in, RField<3>& out);

   /*
   * Get the dimensions of the grid for which this was allocated.
   */
   template <int D>
   void FFT<D>::getMeshDimensions(IntVec<D>& meshDimensions) const
   {
      for (int i = 0; i < D; ++i) {
         meshDimensions[i] = meshDimensions_[i];
      }
   }

}
#include "FFT.tpp"
#endif
