#ifndef PSPC_FFT_H
#define PSPC_FFT_H

/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper for real data.
   *
   * \ingroup Pspc_Field_Module
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
      * Setup grid dimensions, plans and work space.
      *
      * \param rField real data on r-space grid
      * \param kField complex data on k-space grid
      */
      void setup(RField<D>& rField, RFieldDft<D>& kField);

      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * This function computes a scaled Fourier transform, which is obtained
      * by dividing the unscaled transform computed by FFTW by the number of 
      * elements. A scaled copy of the input data is copied to a temporary 
      * real array before being passed to the FFT forward transform function.
      *
      * This function does not overwrite or corrupt the input array.
      * 
      * \param in  array of real values on r-space grid
      * \param out  array of complex values on k-space grid
      */
      void forwardTransform(RField<D> const & in, RFieldDft<D>& out) const;

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * This function computes the same unscaled inverse transform as the
      * FFTW library.
      *
      * NOTE: The inverse transform generally overwrites and corrupts its 
      * input. This is the behavior of the complex-to-real transform of the 
      * underlying FFTW library.  See Sec. 2.3 of the FFTW documentation at
      * https://www.fftw.org/fftw3_doc, One-Dimensional DFTs of Real Data:
      * "...the inverse transform (complex to real) has the side-effect of 
      * overwriting its input array, ..."
      *
      * \param in  array of complex values on k-space grid (overwritten)
      * \param out  array of real values on r-space grid
      */
      void inverseTransform(RFieldDft<D>& in, RField<D>& out) const;

      /**
      * Compute inverse (complex-to-real) Fourier transform without destroying input.
      *
      * \param in  array of complex values on k-space grid (device mem)
      * \param out  array of real values on r-space grid (device mem)
      */
      void inverseTransformSafe(RFieldDft<D> const & kField, RField<D>& rField) const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      IntVec<D> const & meshDimensions() const;

      /** 
      * Has this FFT object been setup?
      */
      bool isSetup() const;

   private:

      /// Private r-space array for performing safe transforms.
      mutable RField<D> rFieldCopy_;

      /// Private k-space array for performing safe transforms.
      mutable RFieldDft<D> kFieldCopy_;

      /// Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      /// Number of points in r-space grid
      int rSize_;

      /// Number of points in k-space grid
      int kSize_;

      /// Pointer to a plan for a forward transform.
      fftw_plan fPlan_;

      /// Pointer to a plan for an inverse transform.
      fftw_plan iPlan_;

      /// Have array dimension and plan been initialized?
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
   * Has this object been setup?
   */
   template <int D>
   inline bool FFT<D>::isSetup() const
   {  return isSetup_; }

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline IntVec<D> const & FFT<D>::meshDimensions() const
   {  return meshDimensions_; }

   #ifndef PSPC_FFT_TPP
   // Suppress implicit instantiation
   extern template class FFT<1>;
   extern template class FFT<2>;
   extern template class FFT<3>;
   #endif

} // namespace Pscf::Pspc
} // namespace Pscf
#endif
