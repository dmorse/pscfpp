#ifndef PSSP_R_FIELD_DFT_H
#define PSSP_R_FIELD_DFT_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Pssp
{

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform of a real field on an FFT mesh.
   *
   * \ingroup Pssp_Field_Module
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
      RFieldDft(const RFieldDft<D>& other);

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
      RFieldDft<D>& operator = (const RFieldDft<D>& other);

      using Field<fftw_complex>::allocate;

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RFieldDft is already allocated.
      *
      * \param meshDimensions vector containing number of grid points in each direction
      */
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Return vector of mesh dimensions by constant reference.
      */
      const IntVec<D>& meshDimensions() const;

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

   };

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RFieldDft<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         if (i < D - 1) {
            size *= meshDimensions[i];
         } else {
            size *= (meshDimensions[i]/2 + 1);
         }
      }
      Field<fftw_complex>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RFieldDft<D>::meshDimensions() const
   {  return meshDimensions_; }


   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RFieldDft<D>::serialize(Archive& ar, const unsigned int version)
   {
      Field<fftw_complex>::serialize(ar, version);
      ar & meshDimensions_;
   }

}
}
#include "RFieldDft.tpp"
#endif
