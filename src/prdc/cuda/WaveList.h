#ifndef PRDC_CUDA_WAVE_LIST_H
#define PRDC_CUDA_WAVE_LIST_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>

#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * Class to calculate and store properties of wavevectors.
   *
   * In particular, minimum images, square norms of wavevectors (kSq), and
   * derivatives of the square norms of wavevectors with respect to the
   * lattice parameters (dKSq) are calculated and stored by this class.
   * 
   * Any time the lattice parameters change the clearUnitCellData() method
   * should be called, which will effectively reset the WaveList object so
   * that the wavevector properties will need to be recalculated before
   * being used.
   * 
   * This object calculates these wavevector properties for a mesh of grid
   * points in k-space. If a calculation only requires real-valued fields, 
   * PSCF uses a reduced-size k-space mesh, as output by cuFFT. However, a 
   * full-sized k-space mesh (the same size as the real-space mesh) is 
   * necessary when dealing with complex-valued fields. The k-space mesh
   * used by a WaveList object is determined by the parameter isRealField,
   * which is assigned in the constructor and cannot later be changed. 
   */
   template <int D>
   class WaveList
   {
   public:

      /**
      * Constructor.
      * 
      * \param isRealField  Will this WaveList be used for real-valued fields?
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

      /**
      * Clear all internal data that depends on lattice parameters.
      *
      * Sets hasKSq_ and hasdKSq_ to false, and sets hasMinImages_ to
      * false if the crystall lattice type has variable angles.
      */
      void clearUnitCellData();

      /**
      * Compute minimum images of wavevectors. (Also calculates kSq.)
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
      * This function uses existing mininum images if they are valid, or
      * recomputes them if necessary.
      */
      void computeKSq();

      /**
      * Compute derivatives of |k|^2 w/ respect to unit cell parameters.
      */
      void computedKSq();

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
      * The array has size kSize, and each element is an IntVec<D>.
      * If isRealField is true, kSize is smaller than the size of the 
      * real-space mesh. Otherwise, it is equal.
      */
      HostDArray< IntVec<D> > const & minImages_h() const;

      /**
      * Get the kSq array on the device by reference.
      * 
      * This method returns an RField in which each element is the square
      * magnitude |k|^2 of a wavevector k in the k-space mesh used for the 
      * DFT. If isRealField is true, this k-space mesh is smaller than the 
      * real-space mesh. Otherwise, it is the same size.
      */
      RField<D> const & kSq() const;

      /**
      * Get derivatives of |k|^2 with respect to lattice parameter i.
      *
      * This method returns an RField in which each element is the 
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
      * each gridpoint. The boolean represents whether the inverse of the
      * wave at the given gridpoint is an implicit wave. Implicit here is
      * used to mean any wave that is outside the bounds of the k-grid.
      * 
      * This method will throw an error if isRealField == false, because
      * there are no implicit inverses in such a case.
      */
      DeviceArray<bool> const & implicitInverse() const;

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
      * Does this WaveList correspond to real-valued fields?
      */ 
      bool isRealField() const
      {  return isRealField_; }

   private:

      /**
      * Array containing minimum images for each wave, stored on device.
      *
      * The array has size kSize * D, where kSize is the number of grid 
      * points in reciprocal space. The array is unwrapped into a linear 
      * array in which the first kSize elements of the array contain the
      * the first coordinate for all minimum image, and so on. If 
      * isRealField is true, kSize is smaller than the size of the 
      * real-space mesh. Otherwise, it is equal.
      */
      DeviceArray<int> minImages_;

      /**
      * Array of IntVec<D> minimum images, stored on the host.
      *
      * Each element of minImageVecs_ contains all D coordinates of the
      * minimum image for a single wavevector, stored on the host as an 
      * IntVec<D>. The array has capacity kSize_. If isRealField is true, 
      * kSize_ is smaller than the size of the real-space mesh. Otherwise, 
      * it is equal.
      */
      mutable
      HostDArray< IntVec<D> > minImages_h_;

      /// Array containing values of kSq_, stored on the device.
      RField<D> kSq_;

      /// Array containing all values of dKSq_, stored on the device.
      DeviceArray<cudaReal> dKSq_;

      /// Array of RFields, where each RField is a slice of the dKSq_ array.
      DArray<RField<D> > dKSqSlices_;

      /// Array indicating whether a given gridpoint has an implicit partner
      DeviceArray<bool> implicitInverse_;

      /**
      * Dimensions of the mesh in reciprocal space.
      * 
      * If isRealField_, the reciprocal-space grid is smaller than the 
      * real-space grid, as output by cuFFT. Otherwise, the two grids
      * are the same size.
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

      /// Has memory been allocated for arrays?
      bool isAllocated_;

      /// Have minimum images been computed (array minImages_) ?
      bool hasMinImages_;

      /// Have minimum image vectors been re-ordered in minImageVecs_) ?
      mutable
      bool hasMinImages_h_;

      /// Has the kSq array been computed?
      bool hasKSq_;

      /// Has the dKSq array been computed?
      bool hasdKSq_;

      /// Will this WaveList be used for real-valued fields?
      bool isRealField_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      /// Access associated UnitCell<D> by reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Access associated Mesh<D> by reference.
      Mesh<D> const & mesh() const
      {  return *meshPtr_; }

   };

   // Get the array of minimum images on the device by reference.
   template <int D>
   inline DeviceArray<int> const & WaveList<D>::minImages_d() const
   {
      UTIL_CHECK(hasMinImages_);
      return minImages_;
   }

   // Get the kSq array on the device by reference.
   template <int D>
   inline RField<D> const & WaveList<D>::kSq() const
   {
      UTIL_CHECK(hasKSq_);
      return kSq_;
   }

   // Get a slice of the dKSq array on the device by reference.
   template <int D>
   inline RField<D> const & WaveList<D>::dKSq(int i) const
   {
      UTIL_CHECK(hasdKSq_);
      return dKSqSlices_[i];
   }

   // Get the implicitInverse array by reference.
   template <int D>
   inline DeviceArray<bool> const & WaveList<D>::implicitInverse() const
   {
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(isRealField_);
      return implicitInverse_;
   }

   #ifndef PRDC_CUDA_WAVE_LIST_TPP
   // Suppress implicit instantiation
   extern template class WaveList<1>;
   extern template class WaveList<2>;
   extern template class WaveList<3>;
   #endif

} // Cuda
} // Prdc
} // Pscf
#endif
