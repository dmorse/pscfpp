#ifndef PRDC_CPU_FFT_H
#define PRDC_CPU_FFT_H

/*
* PSCF Package 
*
* Copyright 2016 - 2024, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/CField.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform wrapper.
   *
   * \ingroup Prdc_Cpu_Module
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
      * Setup grid dimensions, FFT plans and work space.
      *
      * \param meshDimensions Dimensions of real-space grid.
      */
      void setup(IntVec<D> const & meshDimensions);

      // Real Data (Real <-> Complex Transforms)

      /**
      * Compute real-to-complex forward Fourier transform.
      *
      * This function computes a scaled Fourier transform, which is 
      * obtained by dividing the unscaled forward transform computed by 
      * FFTW by the number of grid points elements. A scaled copy of 
      * the input data is copied to a temporary real array before being 
      * passed to the FFTW real-to-complex transform function.
      *
      * This function does not overwrite or corrupt the input array.
      * 
      * \param in  array of real values on r-space grid
      * \param out  array of complex values on k-space grid
      */
      void forwardTransform(RField<D> const & in, RFieldDft<D>& out) const;

      /**
      * Compute complex-to-real to real Fourier transform.
      *
      * This function computes the same unscaled inverse transform as the
      * complex-to-real transforms of the FFTW library.
      *
      * NOTE: The inverse transform generally overwrites and corrupts its 
      * input. This is the behavior of the complex-to-real (c2r) transform 
      * of the underlying FFTW library.  See Sec. 2.3 of the FFTW 3
      * documentation at https://www.fftw.org/fftw3_doc, 
      * One-Dimensional DFTs of Real Data: "...the inverse transform 
      * (complex to real) has the side-effect of overwriting its input 
      * array, ..."
      *
      * \param in  array of complex values on k-space grid (overwritten)
      * \param out  array of real values on r-space grid
      */
      void inverseTransform(RFieldDft<D>& in, RField<D>& out) const;

      /**
      * Compute complex-to-real inverse transform without destroying input.
      *
      * This function makes a copy of the input data and passes the copy
      * to the FFTW complex-to-real transform to avoid overwriting the
      * input data.
      *
      * \param in  array of complex values on k-space grid 
      * \param out  array of real values on r-space grid 
      */
      void inverseTransformSafe(RFieldDft<D> const & in, RField<D>& out) 
      const;

      // Complex Data (Complex <-> Complex Transforms)

      /**
      * Compute complex-to-complex forward Fourier transform.
      *
      * This function computes a scaled forward Fourier transform, which 
      * is obtained by dividing the unscaled forward transform computed 
      * by FFTW by the number of grid points elements. A scaled copy of 
      * the input data is copied to a temporary real array before being 
      * passed to the FFTW real-to-complex transform function.
      *
      * This function does not overwrite or corrupt the input array.
      * 
      * \param in  array of complex values on r-space grid
      * \param out  array of complex values on k-space grid
      */
      void forwardTransform(CField<D> const & in, CField<D>& out) const;

      /**
      * Compute complex-to-complex inverse Fourier transform.
      *
      * This function computes the same unscaled inverse transform as the
      * complex-to-complex transforms of the FFTW library.
      *
      * \param in  array of complex values on k-space grid 
      * \param out  array of complex values on r-space grid
      */
      void inverseTransform(CField<D>& in, CField<D>& out) const;

      // Accessors

      /**
      * Return the dimensions of the grid for which this was setup.
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

      /// Private complex field for performing safe transforms.
      mutable CField<D> cFieldCopy_;

      /// Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      /// Number of points in r-space grid (real or complex data)
      int rSize_;

      /// Number of points in k-space DFT grid for real data
      int kSize_;

      /// Plan for a real-to-complex forward transform.
      fftw_plan rcfPlan_;

      /// Plan for a complex-to-real inverse transform.
      fftw_plan criPlan_;

      /// Plan for a forward complex-to-complex transform.
      fftw_plan ccfPlan_;

      /// Plan for an inverse complex-to-complex transform.
      fftw_plan cciPlan_;

      /// Have array dimensions and plans been initialized?
      bool isSetup_;

      /**
      * Make FFTW plans for transform and inverse transform.
      */
      void makePlans(RField<D>& rField, RFieldDft<D>& kField,
                     CField<D>& cFieldIn, CField<D>& cFieldOut);

   };

   // Declarations of explicit specializations

   template <>
   void FFT<1>::makePlans(RField<1>& rField, RFieldDft<1>& kField,
                          CField<1>& cFieldIn, CField<1>& cFieldOut);

   template <>
   void FFT<2>::makePlans(RField<2>& rField, RFieldDft<2>& kField,
                          CField<2>& cFieldIn, CField<2>& cFieldOut);

   template <>
   void FFT<3>::makePlans(RField<3>& rField, RFieldDft<3>& kField,
                          CField<3>& cFieldIn, CField<3>& cFieldOut);

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

   #ifndef PRDC_CPU_FFT_TPP
   // Suppress implicit instantiation
   extern template class FFT<1>;
   extern template class FFT<2>;
   extern template class FFT<3>;
   #endif

} // namespace Pscf::Prdc::Cpu
} // namespace Pscf::Prdc
} // namespace Pscf
#endif
