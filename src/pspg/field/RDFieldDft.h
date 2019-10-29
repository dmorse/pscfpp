#ifndef PSPG_R_DFIELD_DFT_H
#define PSPG_R_DFIELD_DFT_H

/*
* PSCF++ Package 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DField.h"
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>

namespace Pscf {
namespace Pspg
{

   using namespace Util;
   using namespace Pscf;

   /**
   * Fourier transform of a real field on an FFT mesh.
   *
   * \ingroup Pspg_Field_Module
   */
   template <int D>
   class RDFieldDft : public DField<cufftComplex>
   {

   public:

      /**
      * Default constructor.
      */
      RDFieldDft();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RFieldDft to be copied.
      */
      RDFieldDft(const RDFieldDft<D>& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RDFieldDft();

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
      RDFieldDft<D>& operator = (const RDFieldDft<D>& other);

      using DField<cufftComplex>::allocate;

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
   void RDFieldDft<D>::allocate(const IntVec<D>& meshDimensions)
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
      DField<cufftComplex>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RDFieldDft<D>::meshDimensions() const
   {  return meshDimensions_; }


   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RDFieldDft<D>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent Field capacities");
            }
         }
      }

      if (isAllocated()) {
         cufftComplex* tempData = new cufftComplex[capacity];
         cudaMemcpy(tempData, data_, capacity * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i].x;
            ar & tempData[i].y;
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

}
}
#include "RDFieldDft.tpp"
#endif
