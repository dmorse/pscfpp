#ifndef PRDC_CUDA_R_FIELD_H
#define PRDC_CUDA_R_FIELD_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "types.h"
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Field of real values on a regular mesh, allocated on a GPU device.
   *
   * Type cudaReal is float or double, depending on preprocessor macro.
   *
   * \ingroup Prdc_Cuda_Module 
   */
   template <int D>
   class RField : public DeviceArray<cudaReal>
   {

   public:

      /**
      * Default constructor.
      */
      RField();

      /**
      * Allocating constructor.
      *
      * Allocates memory by calling allocate(meshDimensions) internally.
      *  
      * \param meshDimensions number of grid points in each dimension
      */
      RField(IntVec<D> const & meshDimensions);

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RField to be copied.
      */
      RField(RField<D> const& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RField();

      /**
      * Allocate the underlying C array for data on a regular mesh.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param meshDimensions number of grid points in each dimension
      */
      void allocate(IntVec<D> const & meshDimensions);

      /**
      * Associate this object with a slice of another DeviceArray.
      *
      * \throw Exception if the array is already allocated.
      *
      * \param arr parent array that owns the data
      * \param beginId index in the parent array at which this array starts
      * \param meshDimensions number of grid points in each dimension
      */
      void associate(DeviceArray<cudaReal>& arr, int beginId, 
                     IntVec<D> const & meshDimensions);

      /**
      * Assignment operator, assignment from another RField<D>.
      *
      * Performs a deep copy, by copying all elements of the RHS RField<D>
      * from device memory to device memory, and copying the 
      * meshDimensions.
      *
      * The RHS RField<D> must be allocated on entry. If this LHS object is 
      * not allocated, allocate with the required capacity.  If the LHS and
      * RHS arrays are both allocated, capacity values must be equal.
      * 
      * \param other the RHS RField
      */
      RField<D>& operator = (const RField<D>& other);

      /**
      * Assignment operator, assignment from a HostDArray<cudaReal>.
      *
      * Performs a deep copy, by copying all elements of the RHS RField<D>
      * from host memory to device memory.
      *
      * The RHS HostDArray<cudaReal> and LHS RField<D> must both be 
      * allocated with equal capacity values on entry. 
      * 
      * \param other the RHS HostDArray<cudaReal>
      */
      RField<D>& operator = (const HostDArray<cudaReal>& other);

      /**
      * Return mesh dimensions by constant reference.
      */
      IntVec<D> const & meshDimensions() const;

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
      using DeviceArray<cudaReal>::allocate;

      // Make private to prevent association without mesh dimensions.
      using DeviceArray<cudaReal>::associate;

      // Make private to prevent assignment without mesh dimensions.
      using DeviceArray<cudaReal>::operator=;

   };

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
         HostDArray<cudaReal> tempData(capacity);
         tempData = this; // copy this object's data from device to host
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i];
         }
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
#endif
