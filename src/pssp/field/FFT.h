#ifndef PSSP_FFT_H
#define PSSP_FFT_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp/field/RField.h>
#include <pssp/field/RFieldDft.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Pssp {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Pssp_Field_Module
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
      * Check and setup grid dimensions if necessary.
      *
      * \param rField real data on r-space grid
      * \param kField complex data on k-space grid
      */
      void setup(RField<D>& rField, RFieldDft<D>& kField);

      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-space grid
      * \param out  array of complex values on k-space grid
      */
      void forwardTransform(RField<D>& in, RFieldDft<D>& out);

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid
      * \param out  array of real values on r-space grid
      */
      void inverseTransform(RFieldDft<D>& in, RField<D>& out);

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

   private:

      // Work array for real data.
      RField<D> work_;

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Number of points in r-space grid
      int rSize_;

      // Number of points in k-space grid
      int kSize_;

      // Pointer to a plan for a forward transform.
      fftw_plan fPlan_;

      // Pointer to a plan for an inverse transform.
      fftw_plan iPlan_;

      // Have array dimension and plan been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans for transform and inverse transform.
      */
      void makePlans(RField<D>& rField, RFieldDft<D>& kField);

   };

   // Declarations of explicit specializations

   template <>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField);

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField);

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField);

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline const IntVec<D>& FFT<D>::meshDimensions() const
   {  return meshDimensions_; }

}
}
#include "FFT.tpp"
#endif
