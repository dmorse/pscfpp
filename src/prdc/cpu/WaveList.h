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
   * A WaveList computes and stores minimum images of wavevectors,
   * square norms of wavevectors (kSq), and derivatives of the square
   * norms with respect to the unit cell parameters (dKSq).
   *
   * Any time the lattice parameters change the clearUnitCellData() method
   * should be called. This function sets internal flags that mark some
   * properties as being outdated, indicating that they should recalculated
   * before the next use.
   *
   * A WaveList computes these properties for a mesh of grid points
   * in k-space. If a calculation only involves real-valued fields,
   * PSCF uses a reduced-size k-space mesh, as output by FFTW for the
   * result of a real-to-complex discrete Fourier transform. However, a
   * full-sized k-space mesh (the same size as the real-space mesh) is
   * necessary when dealing with complex-valued fields. The choice of
   * which k-space mesh used by a WaveList object is determined by the
   * bool parameter isRealField that is passed to the constructor, which
   * cannot be changed after construction.
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
      * Sets hasKSq_ and hasdKSq_ to false. Sets hasMinImages_ to
      * false only if the unit cell type has variable angles.
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
      * This function updates the sortedIds and bunchIds arrays.
      */
      void sortWaves();

      ///@}
      /// \name Data Access
      ///@{

      /**
      * Get the array of minimum image vectors by const reference.
      *
      * This function returns an array of kSize elementw in which each 
      * element is an IntVec<D> containing the integer coordinates of 
      * the minimum image of one wavevector in the k-space mesh used 
      * for discrete Fourier transforms. If isRealField is
      * true, this k-space mesh is smaller than the real-space mesh.
      * Otherwise, it is the same size.
      */
      DArray< IntVec<D> > const & minImages() const;

      /**
      * Get the kSq array on the device by const reference.
      *
      * This function returns an array in which each element is the square
      * magnitude |k|^2 of a wavevector k in the k-space mesh used for the
      * DFT. If isRealField is true, this k-space mesh is smaller than the
      * real-space mesh. Otherwise, it is the same size.
      */
      RField<D> const & kSq() const;

      /**
      * Get derivatives of kSq with respect to unit cell parameter i.
      *
      * Each element contains the derivative of the square-wavevector with
      * respect to unit cell parameter i, multiplied by a prefactor. The
      * prefactor is 2.0 for waves that have an implicit inverse and 1.0
      * otherwise. The choice of prefactor is designed to simplify use
      * of the array to compute stress.
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
      * Get all derivatives of kSq with respect to unit cell parameters.
      *
      * Element i of the DArray is the RField<D> that can also be obtained
      * from member function dKSq(int i).
      */
      DArray< RField<D> > const & dKSq() const;

      /**
      * Get the implicitInverse array by reference.
      *
      * This array is defined on a k-grid mesh, with a boolean value for
      * each gridpoint. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave. Implicit here is
      * used to mean any wave that is outside the bounds of the k-grid.
      *
      * This method will throw an error if isRealField == false, because
      * there are no implicit inverses in such a case.
      */
      DArray<bool> const & implicitInverse() const;

      /**
      * Get the sortedIds array by reference.
      *
      * This method will throw an Exception if isSorted == false.
      */
      DArray<int> const & sortedIds() const;

      /**
      * Get the bunchIds array by reference.
      *
      * This method will throw an Exception if isSorted == false.
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
      * Array of indices of bunches to which each wave belongs.
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
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(isSorted_);
      return sortedIds_;
   }

   // Get the bunchIds array by const reference.
   template <int D>
   inline
   DArray<int> const & WaveList<D>::bunchIds() const
   {
      UTIL_CHECK(isAllocated_);
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
