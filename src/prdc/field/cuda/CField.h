#ifndef PRDC_C_FIELD_CU_H
#define PRDC_C_FIELD_CU_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/CUT.h>        // class template argument
#include <pscf/cuda/DeviceArray.h>   // base class template
#include <pscf/cuda/cudaTypes.h>     // base class argument

#include <pscf/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class CField;

   /**
   * Field of complex values on a regular mesh, allocated on a GPU device.
   *
   * \ingroup Prdc_Cuda_Module 
   */
   template <int D>
   class CField<D,CUT>
    : public DeviceArray<cudaComplex>
   {
   
   public:

      /**
      * Default constructor.
      */
      CField();

      /**
      * Allocating constructor.
      *
      * Calls allocate(meshDimension) internally. 
      *
      * \param meshDimensions numbers of grid points in each dimension
      */
      CField(IntVec<D> const & meshDimensions);

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      * \param other the CField to be copied.
      */
      CField(CField const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~CField();

      /**
      * Assignment operator, assignment from another CField<D,CUT>.
      *
      * Performs a deep copy, by copying all elements of the RHS 
      * CField<D,CUT> from device memory to device memory.
      *
      * The RHS CField must be allocated. If this LHS CField is not 
      * allocated on entry, allocate it before copying elements. If 
      * both LHS and RHS objects are allocated on entry, the 
      * capacities must be equal. 
      * 
      * \param other the RHS CField<D,CUT>
      */
      CField<D,CUT>& operator = (CField<D,CUT> const & other);

      /**
      * Assignment operator, assignment from a HostDArray<cudaComplex>.
      *
      * Performs a deep copy, by copying all elements of the RHS 
      * CField<D,CUT> from host memory to device memory.
      *
      * The RHS HostDArray<cudaComplex> and LHS CField<D,CUT> must both be 
      * allocated with equal capacity values on entry. 
      * 
      * \param other the RHS HostDArray<cudaComplex>
      */
      CField<D,CUT>& operator = (HostDArray<cudaComplex> const & other);

      /**
      * Allocate the underlying C array for data on a regular mesh.
      *
      * \throw Exception if the CField is already allocated.
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
      void associate(DeviceArray<cudaComplex>& arr, int beginId, 
                     IntVec<D> const & meshDimensions);

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

      // Vector containing number of grid points in each dimension.
      IntVec<D> meshDimensions_;

      // Make private to prevent allocation without mesh dimensions.
      using DeviceArray<cudaComplex>::allocate;

      // Make private to prevent association without mesh dimensions.
      using DeviceArray<cudaComplex>::associate;

      // Make private to prevent assignment without mesh dimensions.
      using DeviceArray<cudaComplex>::operator=;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& CField<D,CUT>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void CField<D,CUT>::serialize(Archive& ar, const unsigned int version)
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
         HostDArray<cudaComplex> tempData(capacity);
         tempData = this; // copy this object's data from device to host
         for (int i = 0; i < capacity_; ++i) {
            ar & tempData[i].x;
            ar & tempData[i].y;
         }
      }
      ar & meshDimensions_;
   }

   extern template class CField<1,CUT>;
   extern template class CField<2,CUT>;
   extern template class CField<3,CUT>;

}
}
#endif
