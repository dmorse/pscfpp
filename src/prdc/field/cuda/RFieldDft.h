#ifndef PRDC_R_FIELD_DFT_CU_H
#define PRDC_R_FIELD_DFT_CU_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/CUT.h>
#include <pscf/backend/cuda/DeviceArray.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class RFieldDft;

   /**
   * Discrete Fourier Transform (DFT) of a real field, allocated on a GPU.
   *
   * The DFT is stored internally as a C array of cudaComplex elements 
   * located in global GPU memory. All member functions are C++ functions 
   * that can be called from the host CPU. 
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class RFieldDft<D,CUT>
    : public DeviceArray<cudaComplex>
   {

   public:

      /**
      * Default constructor.
      */
      RFieldDft();

      /**
      * Allocating constructor.
      *
      * Calls allocate(meshDimensions) internally.
      *
      * \param meshDimensions dimensions of asssociated r-grid mesh
      */
      RFieldDft(IntVec<D> const & meshDimensions);

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RFieldDft to be copied.
      */
      RFieldDft(RFieldDft<D,CUT> const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RFieldDft();

      /**
      * Assignment operator, assignment from another RFieldDft<D,CUT>.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS Field
      */
      RFieldDft<D,CUT>& operator = (RFieldDft<D,CUT> const & other);

      /**
      * Assignment operator, assignment from a HostDArray<cudaComplex>.
      *
      * Performs a deep copy, by copying all elements of the RHS 
      * RFieldDft<D,CUT> from host memory to device memory.
      *
      * The RHS HostDArray<cudaComplex> and LHS RFieldDft<D,CUT> must 
      * both be allocated and have equal capacity values on entry. 
      * 
      * \param other the RHS HostDArray<cudaComplex>
      */
      RFieldDft<D,CUT>& operator = (HostDArray<cudaComplex> const & other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RFieldDft is already allocated.
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
      * Return vector of real-space mesh dimensions by constant reference.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Return vector of dft (Fourier) grid dimensions by const reference.
      *  
      * The last element of dftDimensions() and meshDimensions() differ by
      * about a factor of two: 
      *    dftDimensions()[D-1] = meshDimensions()/2 + 1.
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

      // Make private to prevent allocation without setting meshDimensions.
      using DeviceArray<cudaComplex>::allocate;

      // Make private to prevent association without setting meshDimensions.
      using DeviceArray<cudaComplex>::associate;

      // Make private to prevent assignment without setting meshDimensions.
      using DeviceArray<cudaComplex>::operator =;

   };

   // Inline and templated member functions
   
   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D> inline 
   const IntVec<D>& RFieldDft<D,CUT>::meshDimensions() const
   {  return meshDimensions_; }

   /*  
   * Return dimensions of dft grid by constant reference. 
   */
   template <int D> inline 
   IntVec<D> const & RFieldDft<D,CUT>::dftDimensions() const
   {  return dftDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RFieldDft<D,CUT>::serialize(Archive& ar, 
		                            const unsigned int version)
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

   extern template class RFieldDft<1,CUT>;
   extern template class RFieldDft<2,CUT>;
   extern template class RFieldDft<3,CUT>;

} // namespace Prdc
} // namespace Pscf
#endif
