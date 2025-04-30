#ifndef PRDC_CUDA_FFT_H
#define PRDC_CUDA_FFT_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RField.h"
#include "CField.h"
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
   * Fourier transform wrapper for real or complex data.
   *
   * This class is a wrapper for plan creation and discrete Fourier 
   * transform (DFT) functions provided by the NVIDIA cufft library, 
   * providing an interface to the field container classes RField<D>, 
   * RField<Dft>, and CField<D> in namespace Pscf::Prdc::Cuda.
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

      // Real <-> Complex transforms
      
      /**
      * Compute forward (real-to-complex) discrete Fourier transform.
      * 
      * The resulting complex array is the output of cuFFT's forward 
      * transform method scaled by a factor of 1/N (where N is the
      * number of grid points). The scaling ensures that a round-trip
      * Fourier transform (forward + inverse) reproduces the original
      * input array.
      * 
      * This transform is inherently "safe", meaning that the input array 
      * will not be modified during the DFT.
      *
      * \param rField  real values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransform(RField<D> const & rField, RFieldDft<D>& kField) 
      const;

      /**
      * Compute inverse (complex-to-real) DFT, overwriting the input.
      * 
      * The resulting real array is the unaltered output of cuFFT's inverse 
      * transform function.
      * 
      * This transform is inherently "unsafe", meaning that the input array
      * will be overwritten during the DFT. Accordingly, the input array
      * kField is not const like it is in all the other transform methods.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransformUnsafe(RFieldDft<D>& kField, RField<D>& rField) 
      const;

      /**
      * Compute inverse (complex-to-real) DFT without overwriting input.
      * 
      * This method makes a copy of the input kField array before performing
      * the DFT, thus preserving the input.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  real values on r-space grid (output, gpu mem)
      */
      void inverseTransformSafe(RFieldDft<D> const & kField, 
                                RField<D>& rField) const;

      // Complex <-> Complex transforms
      
      /**
      * Compute forward (complex-to-complex) discrete Fourier transform.
      * 
      * The resulting complex array is the output of cuFFT's forward 
      * transform method scaled by a factor of 1/N (where N is the
      * number of grid points). The scaling ensures that a round-trip
      * Fourier transform (forward + inverse) reproduces the original
      * input array.
      * 
      * This transform is inherently "safe", meaning that the input array 
      * will not be modified during the DFT.
      * 
      * Note that, in this transform, the k-space grid is the same size
      * as the r-space grid, which differs from the real-to-complex 
      * transform.
      *
      * \param rField  complex values on r-space grid (input, gpu mem)
      * \param kField  complex values on k-space grid (output, gpu mem)
      */
      void forwardTransform(CField<D> const & rField, CField<D>& kField) 
      const;

      /**
      * Compute inverse (complex-to-complex) discrete Fourier transform.
      * 
      * The resulting complex array is the unaltered output of cuFFT's 
      * inverse transform function.
      * 
      * This transform is inherently "safe", meaning that the input array 
      * will not be modified during the DFT.
      * 
      * Note that, in this transform, the k-space grid is the same size
      * as the r-space grid, which differs from the real-to-complex 
      * transform.
      *
      * \param kField  complex values on k-space grid (input, gpu mem)
      * \param rField  complex values on r-space grid (output, gpu mem)
      */
      void inverseTransform(CField<D> const & kField, CField<D>& rField) 
      const;

      /**
      * Return the dimensions of the grid for which this was allocated.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      *  Has this FFT object been setup?  
      */
      bool isSetup() const;

      // Static function

      /**
      * Compute dimensions and size of k-space mesh for DFT of real data.
      * 
      * A corresponding function is not needed for complex-to-complex
      * transforms because the real-space and Fourier-space grids have
      * the same dimensions in this case.
      *
      * \param rMeshDimensions  dimensions of real space grid (real data)
      * \param kMeshDimensions  dimensions of k-space grid (complex data)
      * \param kSize  number of point in k-space grid
      */
      static 
      void computeKMesh(IntVec<D> const & rMeshDimensions,
                        IntVec<D> & kMeshDimensions,
                        int & kSize );

      /**
      * Does a wavevector have an implicit inverse in DFT for real data?
      *
      * The inverse of a wavevector within the DFT mesh used for a DFT of
      * real data is "implicit" if its inverse is not equivalent to any
      * wavevector in this k-mesh.  
      *
      * To compute a Fourier sum using the DFT of real r-space data, for
      * any summand that has inversion symmetry (i.e., the same value
      * for every vector and its inverse), sum over all vectors in the
      * k-space mesh and multiply each element for each wavevector that
      * has an implicit inverse by 2.0, and all others by 1.0. 
      *
      * \param wavevector  integer indices of wavevector (in) 
      * \param meshDimensions  dimensions of real-space mesh (in)
      * \return true iff inverse is implicit (i.e., not in DFT k-mesh)
      */
      static 
      bool hasImplicitInverse(IntVec<D> const & wavevector,
                              IntVec<D> const & meshDimensions);

   private:

      /// Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      /// Private RFieldDft<D> k-space array work space
      mutable RFieldDft<D> kFieldCopy_;

      /// Number of points in r-space grid
      int rSize_;

      /// Number of points in k-space grid for transform of real data
      int kSize_;

      /// Pointer to a plan for a real-to-complex forward transform
      cufftHandle rcfPlan_;

      /// Pointer to a plan for a complex-to-real inverse transform
      cufftHandle criPlan_;

      /// Pointer to a plan for a complex-to-complex transform
      cufftHandle ccPlan_;

      /// Have array dimension and plans been initialized?
      bool isSetup_;

      /**
      * Make cuFFT plans all transform types.
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

   // Inline functions

   /*
   * Return the dimensions of the grid for which this was allocated.
   */
   template <int D>
   inline IntVec<D> const & FFT<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Has this FFT object been setup? 
   */
   template <int D>
   inline bool FFT<D>::isSetup() const
   {  return isSetup_; }

   /*
   * Does this wavevector have an implicit inverse?
   */
   template <int D>
   inline
   bool FFT<D>::hasImplicitInverse(IntVec<D> const & wavevector,
                                   IntVec<D> const & meshDimensions)
   {
      int i = wavevector[D-1];
      int d = meshDimensions[D-1];
      if ((i != 0) && (d - i > d/2 + 1)) {
         return true;
      } else {
         return false;
      }
   }

   #ifndef PRDC_CUDA_FFT_TPP
   // Suppress implicit instantiation
   extern template class FFT<1>;
   extern template class FFT<2>;
   extern template class FFT<3>;
   #endif

} // Cuda
} // Prdc
} // Pscf
#endif
