#ifndef PRDC_CUDA_R_FIELD_H
#define PRDC_CUDA_R_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

#include <cufft.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of real single precision values on an FFT mesh on a device.
   *
   * cudaReal = float or double, depending on preprocessor macro.
   *
   * \ingroup Prdc_Cuda_Module 
   */
   template <int D>
   class RField : public Field<cudaReal>
   {

   public:

      /**
      * Default constructor.
      */
      RField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RField to be copied.
      */
      RField(const RField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RField();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, launch a kernel to swap memory.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      * 
      * \param other the RHS RField
      */
      RField& operator = (const RField& other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param meshDimensions number of grid points in each direction
      */
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Return mesh dimensions by constant reference.
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

      using Field<cudaReal>::allocate;
      using Field<cudaReal>::operator=;

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

   };

   /*
   * Allocate the underlying C array for an FFT grid.
   */
   template <int D>
   void RField<D>::allocate(const IntVec<D>& meshDimensions)
   {
      int size = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(meshDimensions[i] > 0);
         meshDimensions_[i] = meshDimensions[i];
         size *= meshDimensions[i];
      }
      Field<cudaReal>::allocate(size);
   }

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RField<D>::serialize(Archive& ar, const unsigned int version)
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
         double* tempData = new double[capacity];
         cudaMemcpy(tempData, data_, capacity * sizeof(cudaReal), 
                    cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i];
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

   #ifndef PRDC_CUDA_R_FIELD_TPP
   extern template class RField<1>;
   extern template class RField<2>;
   extern template class RField<3>;
   #endif


}
}
}
//#include "RField.tpp"
#endif
