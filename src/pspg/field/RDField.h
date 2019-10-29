#ifndef PSPG_R_DFIELD_H
#define PSPG_R_DFIELD_H

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
   * Field of real single precision values on an FFT mesh on a device.
   *
   * cufftReal = float
   *
   * \ingroup Pspg_Field_Module 
   */
   template <int D>
   class RDField : public DField<cufftReal>
   {

   public:

      /**
      * Default constructor.
      */
      RDField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *uses memcpy! slow!
      *\param other the RField to be copied.
      */
      RDField(const RDField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RDField();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, launch a kernel to swap memory.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      * 
      * uses memcpy! slow!
      * \param other the RHS RField
      */
      RDField& operator = (const RDField& other);

      using DField<cufftReal>::allocate;

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param meshDimensions vector containing number of grid points in each direction
      */
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Return mesh dimensions by constant reference.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      * Serialize a Field to/from an Archive.
      * Temporarily uses a memcpy
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
   void RDField<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      DField<cufftReal>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RDField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RDField<D>::serialize(Archive& ar, const unsigned int version)
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
         float* tempData = new float[capacity];
         cudaMemcpy(tempData, data_, capacity * sizeof(cufftReal), cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i];
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

      



}
}
#include "RDField.tpp"
#endif
