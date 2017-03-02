#ifndef CYLN_FFT_H
#define CYLN_FFT_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/RArray.h>
#include <util/containers/Array.h>
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
      void setup(Array<double>& rField, Array<fftw_complex>& kField);

      /**
      * Compute forward (real-to-complex) Fourier transform.
      *
      * \param in  array of real values on r-space grid
      * \param out  array of complex values on k-space grid
      */
      void forwardTransform(Array<double>& in, Array<fftw_complex>& out);

      /**
      * Compute inverse (complex-to-real) Fourier transform.
      *
      * \param in  array of complex values on k-space grid
      * \param out  array of real values on r-space grid
      */
      void inverseTransform(Array<fftw_complex>& in, Array<double>& out);

   private:

      // Work array for real data.
      Array<double> work_;

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

   };

}
}
#include "FFT.tpp"
#endif
