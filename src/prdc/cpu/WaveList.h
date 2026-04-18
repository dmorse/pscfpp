#ifndef PRDC_CPU_WAVE_LIST_H
#define PRDC_CPU_WAVE_LIST_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/RField.h>           // member
#include <pscf/math/IntVec.h>          // member
#include <util/containers/DArray.h>    // member
#include <util/containers/GArray.h>    // member
#include <util/containers/Pair.h>      // member

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Prdc {
namespace Cpu {

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
   * \ingroup Prdc_Cpu_Module
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
      * Sets hasKSq_ and hasdKSq_ to false. Sets hasMinImages_ to false
      * if the unit cell type has variable angles. Sets isSorted_ to 
      * false if the unit cell has more than one unit cell parameter.
      */
      void clearUnitCellData();

      /**
      * Compute minimum images of wavevectors, and also calculates kSq.
      *
      * This function recomputes the minimum images of all wavevectors if
      * necessary (i.e., if hasMinImages() == false), but does nothing if
      * if minimum images are up to date (if hasMinImages() == true).
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
      * nothing if these values are up to date (if hasKSq() == true).
      * Minimum image values are updated if necessary.
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
      * Get the array of minimum image vectors by const reference.
      *
      * This function returns an array of kSize elements in which each 
      * element is an IntVec<D> containing the integer coordinates of 
      * the minimum image of one wavevector in the k-space mesh used for
      * a discrete Fourier transforms. Array indices correspond to 
      * ranks in this k-space mesh.
      */
      DArray< IntVec<D> > const & minImages() const;

      /**
      * Get the kSq array on the device by const reference.
      *
      * This function returns an array in which each element is the square
      * magnitude |k|^2 of a wavevector in the k-space mesh used for a
      * discrete Fourier transform. Array indices correspond to ranks
      * within this k-space mesh.
      */
      RField<D> const & kSq() const;

      /**
      * Get derivatives of kSq with respect to unit cell parameter i.
      *
      * This function returns an array in which element j contains the
      * derivative of the square-wavevector kSq[j] wtih respect to unit
      * cell parameter number i, multiplied by a weight factor. If the
      * flag isRealField is true, then the weight factor is 2.0 for waves
      * that have an implicit inverse and and 1.0 otherwise. If isReaField
      * is false, then the weight factor is 1.0 for all wavevectors. The 
      * inclusion of a weight factor is designed to simplify use of this 
      * array to compute stress. 
      *
      * Values of the grid array index j correspond to the the rank of 
      * the k-space mesh used for discrete Fourier transforms.
      *
      * \param i index of lattice parameter
      */
      RField<D> const & dKSq(int i) const;

      /**
      * Get all derivatives of kSq with respect to unit cell parameters.
      *
      * Element i of the DArray is the RField<D> that can also be obtained
      * from member function dKSq(int i). See documentation of that
      * function.
      */
      DArray< RField<D> > const & dKSq() const;

      /**
      * Get the implicitInverse array by reference.
      *
      * This array is defined on a k-grid mesh, with a boolean value for
      * each gridpoint. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave in the k-space
      * mesh used for a discrete Fourier transform of a real function. 
      * The inverse is implicit if it is outside the bounds of this
      * truncated k-space mesh. Array indices correspond to ranks within
      * this k-space mesh.
      *
      * This function throws an Exception if isRealField == false, because
      * there are no implicit inverses in such a case.
      */
      DArray<bool> const & implicitInverse() const;

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
      * that give upper and lower bounds for array indices in the 
      * sortedIds array of a contiguous slice (a "bunch") of sorted
      * waves for which the wavevector have equal vector magitudes. The 
      * first integer in such a pair is the array index of the first 
      * wavevector in such a bunch, and the second is one greater than 
      * the index of the last wavevector in the bunch.
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
      * Return the dimensions of the k-grid mesh.
      * 
      * If isRealField() == true, the reciprocal-space grid is smaller 
      * than the real-space grid. Otherwise, the two grids are identical.
      */
      IntVec<D> const & kMeshDimensions() const
      {  return kMeshDimensions_; }

      /**
      * Return the number of points in the k-grid mesh.
      *
      * If isRealField() == true, kSize is approximately half the size
      * of the real-space grid.  Otherwise, the two grids are identical.
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
      * Are minimum images up to date ?
      */
      bool hasMinImages() const
      {  return hasMinImages_; }

      /**
      * Are values of kSq up-to-date ?
      */
      bool hasKSq() const
      {  return hasKSq_; }

      /**
      * Are values of dKSq up-to-date?
      */
      bool hasdKSq() const
      {  return hasdKSq_; }

      /**
      * Are waves sorted ?
      */
      bool isSorted() const
      {  return isSorted_; }

      /**
      * Does this WaveList correspond to real-valued fields?
      */
      bool isRealField() const
      {  return isRealField_; }

      ///@}

   private:

      /*
      * Array indices for arrays minImages_, kSq_, dKSq_ implicitInverse_,
      * and bunchIds_ correspond to ranks with the mesh with dimensions
      * given by kMeshDimensions_.
      */

      /**
      * Array of minimum images for each wave, indexed by wave rank.
      */
      DArray< IntVec<D> > minImages_;

      /**
      * Array of square-magnitude values for wavevectors.
      */
      RField<D> kSq_;

      /**
      * Derivatives of kSq_ with respect to lattice parameters.
      *
      * Element kSq_[i][j] is the derivative of kSq_[j] with respect to
      * lattice parameter i. 
      */
      DArray< RField<D> > dKSq_;

      /**
      * Array indicating whether a given gridpoint has an implicit partner.
      * 
      * This array is allocated and used only if isRealField == true.
      */
      DArray<bool> implicitInverse_;

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
      * If isRealField_, the reciprocal-space grid is smaller 
      * than the real-space grid, as output by FFTW. Otherwise, the two 
      * grids are identical.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of grid points in reciprocal space.
      *
      * If isRealField_, the reciprocal-space grid is smaller than the
      * real-space grid, as output by cuFFT. Otherwise, the two grids
      * are the same size.
      */
      int kSize_;

      /**
      * Number of distinct wavevector magnitudes.
      */
      int nBunch_;

      /// Has memory been allocated for arrays?
      bool isAllocated_;

      /// Have minimum images been computed?
      bool hasMinImages_;

      /// Has the kSq array been computed?
      bool hasKSq_;

      /// Has the dKSq array been computed?
      bool hasdKSq_;

      /// Have the waves been sorted by magnitude?
      bool isSorted_;

      /// Will this WaveList be used for real-valued fields?
      bool isRealField_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      /// Access associated UnitCell<D> by const reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Access associated Mesh<D> by const reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

   };

   // Get the array of minimum images on the device by reference.
   template <int D>
   inline
   DArray< IntVec<D> > const & WaveList<D>::minImages() const
   {
      UTIL_CHECK(hasMinImages_);
      return minImages_;
   }

   // Get the kSq array on the device by reference.
   template <int D>
   inline
   RField<D> const & WaveList<D>::kSq() const
   {
      UTIL_CHECK(hasKSq_);
      return kSq_;
   }

   // Get dKSq for unit cell parameter array i.
   template <int D>
   inline
   RField<D> const & WaveList<D>::dKSq(int i) const
   {
      UTIL_CHECK(hasdKSq_);
      return dKSq_[i];
   }

   // Get entire dKSq container by const reference.
   template <int D>
   inline
   DArray< RField<D> > const & WaveList<D>::dKSq() const
   {
      UTIL_CHECK(hasdKSq_);
      return dKSq_;
   }

   // Get the implicitInverse array by const reference.
   template <int D>
   inline
   DArray<bool> const & WaveList<D>::implicitInverse() const
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
   template <int D>
   inline
   GArray< Pair<int> > const & WaveList<D>::sortedBunches() const
   {
      UTIL_CHECK(isSorted_);
      return sortedBunches_;
   }

   // Get the bunchIds array by const reference.
   template <int D>
   inline
   DArray<int> const & WaveList<D>::bunchIds() const
   {
      UTIL_CHECK(isSorted_);
      return bunchIds_;
   }

   // Explicit instantiation declarations
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;

} // Cpu
} // Prdc
} // Pscf
#endif
