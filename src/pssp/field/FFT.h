#ifndef PSSP_FFT_H
#define PSSP_FFT_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp/field/RMeshField.h>
#include <pssp/field/KMeshField.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pssp
{

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of real double precision values on an FFT mesh.
   */
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
      template <int D>
      void forwardTransform(RMeshField& in, KMeshField& out);

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid
      * \param out  array of real values on r-space grid
      * \param meshDimensions number of grid points in each direction.
      */
      template <int D>
      void inverseTransform(KMeshField& in, RMeshField& out);

      /**
      * Return the dimension of space.  
      */
      int spaceDimension() const;

      /**
      * Get the dimensions of the grid for which this was allocated.
      *
      * \throw Exception if dimensions of space do not match.
      *
      * \param dimensions number of grid points in each direction.
      */
      template <int D>
      void getMeshDimensions(IntVec<D>& meshDimensions) const;

   private:

      // Vector containing number of grid points in each direction.
      IntVec<3> meshDimensions_;

      // Dimension of space (1, 2, or 3)
      int spaceDimension_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid
      int kSize_;

      // Pointer to a plan for a forward transform.
      fftw_plan fPlan_;

      // Pointer to a plan for a forward transform.
      fftw_plan rPlan_;

      /**
      * Check and setup grid dimensions if necessary.
      *
      * \param dimensions number of grid points in each direction.
      */
      template <int D>
      void setup(const RMeshField& rField, const KMeshField& kField);

      /**
      * Set new grid dimensions.
      *
      * \param dimensions number of grid points in each direction.
      */
      template <int D>
      void setDimensions(const IntVec<D>& meshDimensions);

   };

   // Declarations of explicit specializations

   template <>
   void FFT::forwardTransform<1>(RMeshField& in, KMeshField& out);

   template <>
   void FFT::forwardTransform<2>(RMeshField& in, KMeshField& out);

   template <>
   void FFT::forwardTransform<3>(RMeshField& in, KMeshField& out);

   template <>
   void FFT::inverseTransform<1>(KMeshField& in, RMeshField& out);

   template <>
   void FFT::inverseTransform<2>(KMeshField& in, RMeshField& out);

   template <>
   void FFT::inverseTransform<3>(KMeshField& in, RMeshField& out);

   inline int FFT::spaceDimension() const
   {  return spaceDimension_;}

   /*
   * Get the dimensions of the grid for which this was allocated.
   */
   template <int D>
   void FFT::getMeshDimensions(IntVec<D>& meshDimensions) const
   {
      if (D != spaceDimension_) {
         UTIL_THROW("Argument with wrong number of spatial dimensions");
      } else {
         for (int i = 0; i < D; ++i) {
            meshDimensions[i] = meshDimensions_[i];
         }
      }
   }

   /*
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFT::setup(const RMeshField& rField, const KMeshField& kField)
   {
      IntVec<D> rDimensions;
      IntVec<D> kDimensions;
      rField.getMeshDimensions(rDimensions);
      kField.getMeshDimensions(kDimensions);
      UTIL_CHECK(rDimensions == kDimensions);
      if (spaceDimension_ == 0) {
         setDimensions(rDimensions);
      } else {
         for (int i = 0; i < D; ++i) {
            if (rDimensions[i] != meshDimensions_[i]) {
               UTIL_THROW("Mismatched dimensions");
            }
         }
      }
      UTIL_CHECK(rField.capacity() == rSize_);
      UTIL_CHECK(kField.capacity() == kSize_);
   }

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void FFT::setDimensions(const IntVec<D>& meshDimensions)
   {
      // Preconditions
      UTIL_CHECK(D > 0);
      UTIL_CHECK(D < 4);

      // Initialize mesh dimensions to 1
      for (int i = 0; i < 3; ++i) {
         meshDimensions_[i] = 1;
      }

      int rSize_ = 1;
      int kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         rSize_ *= meshDimensions[i];
         if (i < D - 1) {
            kSize_ *= meshDimensions[i];
         } else {
            kSize_ *= (meshDimensions[i]/2 + 1);
         }
      }
      spaceDimension_ = D;
   }

}
#endif
