#ifndef PRDC_R_FIELD_CU_H
#define PRDC_R_FIELD_CU_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/CUT.h>         // class template argument
#include <pscf/cuda/DeviceArray.h>    // base class template
#include <pscf/cuda/cudaTypes.h>      // base class argument
#include <pscf/math/IntVec.h>         // class member

#include <pscf/cuda/HostDArray.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declaration of primary template
   template <int D, class T> class RField;

   /**
   * Field of real values on a regular mesh, allocated on a GPU device.
   *
   * Type cudaReal is float or double, depending on preprocessor macro.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class RField<D,CUT>
    : public DeviceArray<cudaReal>
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
      RField(RField<D,CUT> const& other);

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
      * Assignment operator, assignment from another RField.
      *
      * Performs a deep copy, by copying all elements of the RHS RField
      * from device memory to device memory, and copying the
      * meshDimensions.
      *
      * The RHS RField must be allocated on entry. If this LHS object is
      * not allocated, allocate with the required capacity.  If the LHS and
      * RHS arrays are both allocated, capacity values must be equal.
      *
      * \param other the RHS RField
      */
      RField<D,CUT>& 
      operator = (RField<D,CUT> const & other);

      /**
      * Assignment operator, assignment from a HostDArray<cudaReal>.
      *
      * Performs a deep copy, by copying all elements of the RHS RField
      * from host memory to device memory.
      *
      * The RHS HostDArray<cudaReal> and LHS RField must both be
      * allocated with equal capacity values on entry.
      *
      * \param other the RHS HostDArray<cudaReal>
      */
      RField<D,CUT>& operator = (const HostDArray<cudaReal>& other);

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
   template <int D> inline 
   IntVec<D> const & 
   RField<D,CUT>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize an RField to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RField<D,CUT>::serialize(Archive& ar, const unsigned int version)
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

   // Explicit instantiation declarations
   extern template class RField<1,CUT>;
   extern template class RField<2,CUT>;
   extern template class RField<3,CUT>;

} // namespace Prdc
} // namespace Pscf
#endif
