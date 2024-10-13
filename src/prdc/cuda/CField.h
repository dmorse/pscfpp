#ifndef PRDC_CUDA_C_FIELD_H
#define PRDC_CUDA_C_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/Field.h>
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   template <typename Data> class HostField;

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of complex values on an regular mesh, allocated on a GPU device.
   *
   * \ingroup Prdc_Cuda_Module 
   */
   template <int D>
   class CField : public Field<cudaComplex>
   {

   public:

      /**
      * Default constructor.
      */
      CField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the CField to be copied.
      */
      CField(const CField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~CField();

      /**
      * Assignment operator, assignment from another CField<D>.
      *
      * Performs a deep copy, by copying all elements of the RHS CField<D>
      * from device memory to device memory.
      *
      * The RHS CField must be allocated. If this LHS CField is not allocated
      * on entry, allocate it before copying elements. If both LHS and RHS
      * objects are allocated on entry, the capacities must be equal. 
      * 
      * \param other the RHS CField<D>
      */
      CField<D>& operator = (const CField<D>& other);

      /**
      * Assignment operator, assignment from a HostField<cudaComplex>.
      *
      * Performs a deep copy, by copying all elements of the RHS CField<D>
      * from host memory to device memory.
      *
      * The RHS HostField<cudaComplex> and LHS CField<D> must both be allocated
      * with equal capacity values on entry. 
      * 
      * \param other the RHS HostField<cudaComplex>
      */
      CField<D>& operator = (const HostField<cudaComplex>& other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the CField is already allocated.
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

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      // Make private to prevent allocation without mesh dimensions.
      using Field<cudaComplex>::allocate;

      // Make private to prevent assignment without mesh dimensions
      using Field<cudaComplex>::operator=;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& CField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void CField<D>::serialize(Archive& ar, const unsigned int version)
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
         cudaMemcpy(tempData, data_, capacity * sizeof(cudaComplex), 
                    cudaMemcpyDeviceToHost);
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i];
         }
         delete[] tempData;
      }
      ar & meshDimensions_;
   }

   #ifndef PRDC_CUDA_C_FIELD_TPP
   extern template class CField<1>;
   extern template class CField<2>;
   extern template class CField<3>;
   #endif


}
}
}
//#include "CField.tpp"
#endif
