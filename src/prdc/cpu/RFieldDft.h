#ifndef PRDC_CPU_R_FIELD_DFT_H
#define PRDC_CPU_R_FIELD_DFT_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/Field.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform of a real field on an FFT mesh.
   *
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class RFieldDft : public Field<fftw_complex>
   {

   public:

      /**
      * Default constructor.
      */
      RFieldDft();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RFieldDft to be copied.
      */
      RFieldDft(RFieldDft<D> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RFieldDft();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS Field
      */
      RFieldDft<D>& operator = (RFieldDft<D> const & other);

      using Field<fftw_complex>::allocate;

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RFieldDft is already allocated.
      *
      * \param meshDimensions vector of grid points in each direction
      */
      void allocate(IntVec<D> const & meshDimensions);

      /**
      * Return vector of spatial mesh dimensions by constant reference.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Return vector of dft (Fourier) grid dimensions by constant reference.
      *  
      * The last element of dftDimensions() and meshDimensions() differ by
      * about a factor of two: dftDimension()[D-1] = meshDimensions()/2 + 1.
      * For D > 1, other elements are equal. 
      */
      IntVec<D> const & dftDimensions() const;

      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Vector containing dimensions of dft (Fourier) grid.
      IntVec<D> dftDimensions_;

   };

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RFieldDft<D>::allocate(IntVec<D> const & meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            dftDimensions_[i] = meshDimensions[i];
            size *= meshDimensions[i];
         } else {
            dftDimensions_[i] = (meshDimensions[i]/2 + 1);
            size *= dftDimensions_[i];
         }
      }
      Field<fftw_complex>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline IntVec<D> const & RFieldDft<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Return dimensions of dft grid by constant reference.
   */
   template <int D>
   inline IntVec<D> const & RFieldDft<D>::dftDimensions() const
   {  return dftDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RFieldDft<D>::serialize(Archive& ar, const unsigned int version)
   {
      Field<fftw_complex>::serialize(ar, version);
      ar & meshDimensions_;
      ar & dftDimensions_;
   }

   #ifndef PRDC_R_FIELD_DFT_TPP
   extern template class RFieldDft<1>;
   extern template class RFieldDft<2>;
   extern template class RFieldDft<3>;
   #endif

}
}
}
#endif
