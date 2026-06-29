#ifndef PRDC_CUDA_WAVE_LIST_H
#define PRDC_CUDA_WAVE_LIST_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/cuda/RField.h>        // member
#include <pscf/cuda/DeviceArray.h>   // member
#include <pscf/cuda/HostDArray.h>    // member
#include <pscf/math/IntVec.h>        // member
#include <util/containers/DArray.h>  // member
#include <util/containers/GArray.h>  // member
#include <util/containers/Pair.h>    // member

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Class to compute and store properties associated with wavevectors.
   *
   * A WaveList computes and stores minimum images of wavevectors, square
   * norms of wavevectors (kSq), and derivatives of the square norms with
   * respect to the unit cell parameters (dKSq). A WaveList can also be 
   * used to sort wavevectors in order of increasing wavevector norm.
   *
   * A WaveList computes these properties for a mesh of grid points in
   * k-space. Different dimensions are used for this mesh depending on the
   * value of the parameter isRealField that is passed to constructor. This
   * parameter specifies whether this class will use a k-space mesh designed
   * for real or complex fields throughout the lifetime of this WaveList 
   * object. If isRealField == true, this class uses the k-space mesh used 
   * by FFTW for the result of a real-to-complex discrete Fourier transform, 
   * which is contains slightly more than half the number of grid points as 
   * the corresponding real-space mesh. If isRealField == true, this class 
   * uses a k-space mesh with dimensions that are the same as those of the 
   * associated real-space mesh. The dimensions of this k-space mesh are
   * returned as an IntVec<D> by member function kMeshDimensions().
   *
   * In most of the arrays constructed by this class, each element 
   * corresponds to a wavevector in the associated k-space mesh, and the
   * index of each element correspond to the rank of the associated 
   * wavevector within that mesh. This applies to the arrays returned by
   * the functions minImages(), kSq(), dkSq(int i), and bunchIds().
   *
   * Any time the lattice parameters change, the clearUnitCellData() method
   * should be called. This function sets internal flags that mark some
   * properties as being outdated, indicating that they are outdated and 
   * thus must be recalculated before the next use.
   *
   * \ingroup Prdc_Cuda_Module
   */
   template <int D>
   class WaveList
   {
   public:

      /// \name Construction, Destruction and Initialization
      ///@{ 

      /**
      * Constructor.
      *
      * \param isRealField  Will this object be used for real-valued fields?
      */
      WaveList(bool isRealField = true);

      /**
      * Destructor
      */
      ~WaveList();

      /**
      * Allocate memory and set association with a Mesh and UnitCell object.
      *
      * \param m  spatial discretization mesh (input)
      * \param c  crystallographic unit cell (input)
      */
      void allocate(Mesh<D> const & m, UnitCell<D> const & c);

      ///@}
      /// \name Computation
      ///@{ 

      /**
      * Clear all internal data that depends on lattice parameters.
      *
      * Sets hasKSq_ and hasdKSq_ to false, and sets hasMinImages_ to
      * false only if the unit cell type has variable angles.
      */
      void clearUnitCellData();

      /**
      * Compute minimum images of wavevectors, and also calculate kSq.
      *
      * This function recomputes the minimum images of all wavevectors if
      * necessary (i.e., if hasMinImages() == false), but does nothing if
      * the minimum images are up to date (if hasMinImages() == true).
      *
      * The minimum images may change if a lattice angle in the unit cell
      * is changed, so this method should be called whenever such changes
      * occur.
      *
      * In the process of computing the minimum images, the square norm
      * |k|^2 for all wavevectors is also calculated and stored, so it
      * is not necessary to call computeKSq after calling this method.
      * computeKSq is provided to allow calculation of kSq without
      * recalculating minimum images.
      */
      void computeMinimumImages();

      /**
      * Compute square norm |k|^2 for all wavevectors.
      *
      * This function recomputes values of the square norm for all 
      * wavevectors if necessary (i.e., if hasKSq() == false), and does
      * nothing if these values are up to date (if hasKSq() ==  true).
      * Mininum image values are updated if necessary. 
      */
      void computeKSq();

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      *
      * This function computes values of the derivatives of wavevector
      * square norms with respect to unit cell parameters if necessary
      * (i.e., if hasdKSq() == false), and does nothing if these values
      * are up to date. Minimum images are updated if necessary.
      */
      void computedKSq();

      /**
      * Sort waves in order of ascending wavevector norm.
      *
      * This function computes the sortedIds, sortedBunches, and bunchIds 
      * arrays. Values of kSq are updated before sorting waves if hasKSq()
      * is false on entry.
      */
      void sortWaves();

      ///@}
      /// \name Data Access
      ///@{ 

      /**
      * Get the array of minimum images on the device by reference.
      *
      * The array has size kSize * D, where kSize is the number of grid
      * points in the FFT k-space mesh. The array is unwrapped into a
      * linear array in an index-by-index manner, in which the first kSize
      * elements of the array contain the first index of each minimum
      * image, and so on. If isRealField is true, kSize is smaller than
      * the size of the real-space mesh. Otherwise, it is equal.
      */
      DeviceArray<int> const & minImages_d() const;

      /**
      * Get minimum images as IntVec<D> objects on the host.
      *
      * This function returns an array of kSize elements in which each
      * element is an IntVec<D> containing the integer coordinates of 
      * the minimum image of one wavevector in the k-space mesh used for
      * discrete Fourier transforms.
      */
      HostDArray< IntVec<D> > const & minImages_h() const;

      /**
      * Get the kSq array on the device by reference.
      *
      * This method returns an RField<D> in which each element is the 
      * square magnitude |k|^2 of a wavevector k in the k-space mesh used 
      * for a DFT. If isRealField is true, this k-space mesh is smaller 
      * than the real-space mesh. Otherwise, it is the same size.
      */
      RField<D> const & kSq() const;

      /**
      * Get derivatives of |k|^2 with respect to lattice parameter i.
      *
      * This method returns an RField<D> in which each element is the
      * derivative of the square-wavevector with respect to unit cell
      * parameter i, multiplied by a prefactor. The prefactor is 2.0 for
      * waves that have an implicit inverse and 1.0 otherwise. The choice
      * of prefactor is designed to simplify use of the array to compute
      * stress.
      *
      * Each element corresponds to one wavevector k in the k-space mesh
      * used for the DFT. If isRealField is true, this k-space mesh is
      * smaller than the real-space mesh. Otherwise, it is the same size.
      * In the latter case, there are no implicit waves, so the prefactor
      * is always 1.0.
      *
      * \param i index of lattice parameter
      */
      RField<D> const & dKSq(int i) const;

      /**
      * Get the implicitInverse array by reference.
      *
      * This array is defined on a k-grid mesh, with a boolean value for
      * each wavevector. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave. Implicit here is
      * used to mean any wave that is outside the bounds of the k-grid.
      *
      * This method will throw an error if isRealField == false, because
      * there are no implicit inverses in such a case.
      */
      DeviceArray<bool> const & implicitInverse() const;

      /**
      * Get the sortedIds array by reference.
      *
      * The array returned by this function contains k-grid rank indices
      * of wavevectors sorted in order of increasing wavevector magnitude.
      *
      * This method throws an Exception if isSorted == false.
      */
      DArray<int> const & sortedIds() const;

      /**
      * Get the sortedBunches array by reference.
      *
      * Each element in this array contains a Pair<int> of two integers 
      * that give upper and lower index bounds of a contiguous slice
      * (a "bunch") that contains ids of wavevectors of equal magnitude.
      * The first value in such a pair is the array index in sortIds_ of 
      * the first element of such bunch, and the second is one greater 
      * than the index of the last element in that bunch.
      *
      * This method throws an Exception if isSorted == false.
      */
      GArray< Pair<int> > const & sortedBunches() const;

      /**
      * Get the bunchIds array by reference.
      *
      * This function returns an array, indexed by wave rank, in which
      * the value is the integer identifier of a bunch to which the wave
      * belongs. Each such value corresponds to an array index of the
      * corresponding bunch in the sortedBunches array.
      *
      * This method throws an Exception if isSorted == false.
      */
      DArray<int> const & bunchIds() const;

      /**
      * Return the dimensions of the k-space mesh.
      *
      * If isRealField() == true, the reciprocal-space grid is smaller
      * than the real-space grid. Otherwise, the two grids are identical.
      */
      IntVec<D> const & kMeshDimensions() const
      {  return kMeshDimensions_; }

      /**
      * Return the number of points in the k-space mesh.
      *
      * If isRealField() == true, kSize is approximately half the size
      * of the real-space grid.  Otherwise, the two grids are identical.
      * 
      */
      int kSize() const
      {  return kSize_; }

      /**
      * Return the number of bunches of sorted waves.
      */
      int nBunch() const
      {
         UTIL_CHECK(isSorted_);  
         return nBunch_; 
      }

      ///@}
      /// \name Boolean Queries
      ///@{ 

      /**
      * Has memory been allocated for arrays?
      */
      bool isAllocated() const
      {  return isAllocated_; }

      /**
      * Have minimum images been computed?
      */
      bool hasMinImages() const
      {  return hasMinImages_; }

      /**
      * Has the kSq array been computed?
      */
      bool hasKSq() const
      {  return hasKSq_; }

      /**
      * Has the dKSq array been computed?
      */
      bool hasdKSq() const
      {  return hasdKSq_; }

      /**
      * Have the waves been sorted by magnitude ?
      */
      bool isSorted() const
      {  return isSorted_; }

      /**
      * Is this WaveList set up for use with real-valued fields?
      */
      bool isRealField() const
      {  return isRealField_; }

      ///@}

   private:

      // Private member variables

      /**
      * Array containing minimum images for all waves, stored on device.
      *
      * The array has size kSize_ * D, where kSize_ is the number of grid
      * points in reciprocal space. The array is unwrapped into a linear
      * array in which the first kSize_ elements of the array contain
      * the the first coordinate for all minimum image, and so on. If
      * isRealField_ is true, kSize_ is smaller than the size of the
      * real-space mesh. Otherwise, kSize_ is equal to the size of the
      * real-space mesh.
      */
      DeviceArray<int> minImages_;

      /**
      * Array of IntVec<D> minimum images for waves, stored on the host.
      *
      * Each element of minImageVecs_ contains all D coordinates of the
      * minimum image for a single wavevector, stored on the host as an
      * IntVec<D>. The array has capacity kSize_. If isRealField is true,
      * kSize_ is smaller than the size of the real-space mesh. Otherwise,
      * kSize_ is equal to the size of the real space mesh.
      */
      mutable
      HostDArray< IntVec<D> > minImages_h_;

      /**
      * Array containing values of kSq_, stored on the device.
      *
      * The mesh dimensions are those of the reciprocal space mesh.
      */
      RField<D> kSq_;

      /**
      * Array containing all values of dKSq_, stored on the device.
      *
      * The dimensions are kSize_ * nParam, where nParam is the number
      * of unit cell parameters.
      */
      DeviceArray<cudaReal> dKSq_;

      /**
      * Array of RFields, where each RField is a slice of the dKSq_ array.
      *
      * The number of elements is equal to nParam, the number of unit cell
      * parameters.  Element dKSqSlices_[i] is an  RField<D> element that
      * is associated with a slice of the larger dKSq_ device array, and
      * that contains derivatives of square wavevectors with respect to
      * unit cell parameter number i.
      *
      * The dKSqSlices_ container should appear after dKSq_ in the
      * declaration of class members in order to gurantee that the
      * elements of dkSqSlices_ will be destroyed before the dKSq_
      * container that owns the data.
      */
      DArray< RField<D> > dKSqSlices_;

      /**
      * Array indicating whether a given gridpoint has an implicit partner.
      *
      * This array is only allocated and used if isRealField_ is true,
      * in which case it has size kSize_.
      */
      DeviceArray<bool> implicitInverse_;

      /**
      * Wavevector ranks, sorted in ascending wavevector magnitude.
      *
      * For i > j, kSq_[sortedIds_[i]] >= kSq_[sortedIds_[j]]. 
      */
      DArray<int> sortedIds_;

      /**
      * Bounds of "bunches" of wavevectors of equal norm in sortedIds_.
      *
      * Each element contains a pair of integers that give bounds of a
      * slice of the array sortedIds_ that contains ids for wavevectors
      * of equal magnitude. The first value in each pair is the index
      * within sortedIds_ of the first element in such a "bunch" and 
      * the second is one greater than the index of the last element 
      * in that bunch.
      */
      GArray< Pair<int> > sortedBunches_;

      /**
      * Array containing ids of bunches to which each wave belongs.
      *
      * The value of element j is the integer identifier of a "bunch"
      * to which wave j belongs. Each such value corresponds to the
      * array index of that bunch in the sortedBunches_ array.
      *
      * If bunchIds_[i] == bunchIds_[j], then kSq_[i] == kSq_[j].
      */
      DArray<int> bunchIds_;

      /**
      * Dimensions of the mesh in reciprocal space.
      *
      * If isRealField_ is true, kMeshDimensions_ is equal to the vector
      * of dimensions of the reciprocal space grid used by a RFieldDft<D>
      * container to store the discrete Fourier transform of a real field.
      * One dimension of this mesh is approximately half the corresponding
      * dimension of the associated real space grid.  If isRealField_
      * is false, indicating application to a complex field, then
      * kMeshDimensions_ is equal to the vector of dimensions of the real
      * space grid.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of grid points in reciprocal space.
      *
      * The integer kSize_ is the number of elements in the reciprocal
      * space grid, given by the product of elements of kMeshDimensions_.
      * If isRealField_, kSize_ is smaller than the size of the real
      * space mesh, by approximately a factor of 2 for large meshes.
      */
      int kSize_;

      /**
      * Number of distinct wavevector magnitudes.
      */
      int nBunch_;

      /// Has memory been allocated for private member arrays?
      bool isAllocated_;

      /// Do valid minimum images exist (array minImages_) ?
      bool hasMinImages_;

      /// Have minimum image vectors been re-ordered in minImageVecs_ ?
      mutable
      bool hasMinImages_h_;

      /// Do valid values of kSq_ array exist ?
      bool hasKSq_;

      /// Do valid values of dKSq_ exist ?
      bool hasdKSq_;

      /// Have the waves been sorted by magnitude?
      bool isSorted_;

      /// Will this WaveList be used for real-valued fields?
      bool isRealField_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      // Private member functions

      /// Access associated UnitCell<D> by reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Access associated Mesh<D> by reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

   };

   // Get the array of minimum images on the device by reference.
   template <int D> inline
   DeviceArray<int> const & WaveList<D>::minImages_d() const
   {
      UTIL_CHECK(hasMinImages_);
      return minImages_;
   }

   // Get the kSq array on the device by reference.
   template <int D> inline
   RField<D> const & WaveList<D>::kSq() const
   {
      UTIL_CHECK(hasKSq_);
      return kSq_;
   }

   // Get a slice of the dKSq array on the device by reference.
   template <int D> inline 
   RField<D> const & WaveList<D>::dKSq(int i) const
   {
      UTIL_CHECK(hasdKSq_);
      return dKSqSlices_[i];
   }

   // Get the implicitInverse array by reference.
   template <int D> inline
   DeviceArray<bool> const & WaveList<D>::implicitInverse() const
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(isRealField_);
      return implicitInverse_;
   }

   // Get the sortedIds array by const reference.
   template <int D>
   inline
   DArray<int> const & WaveList<D>::sortedIds() const
   {
      UTIL_CHECK(isSorted_);
      return sortedIds_;
   }

   // Get the sortedBunches array by const reference.
   template <int D> inline
   GArray< Pair<int> > const & WaveList<D>::sortedBunches() const
   {
      UTIL_CHECK(isSorted_);
      return sortedBunches_;
   }

   // Get the bunchIds array by const reference.
   template <int D> inline
   DArray<int> const & WaveList<D>::bunchIds() const
   {
      UTIL_CHECK(isSorted_);
      return bunchIds_;
   }

   // Explicit instantiation declarations
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif
