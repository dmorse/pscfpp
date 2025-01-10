#ifndef RPG_BLOCK_H
#define RPG_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"                  // base class argument
#include <prdc/cuda/RField.h>            // member
#include <prdc/cuda/RFieldDft.h>         // member
#include <prdc/cuda/FFT.h>               // member
#include <prdc/cuda/FFTBatched.h>        // member
#include <rpg/solvers/WaveList.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/solvers/BlockTmpl.h>      // base class template
#include <util/containers/FSArray.h>

namespace Pscf {

   template <int D> class Mesh;

namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   // Wrappers for the CUDA kernels used internally
   
   /**
   * Element-wise calculation of a = real(b * conj(c) * d), kernel wrapper
   * 
   * \param a  output array (real)
   * \param b  input array 1 (complex)
   * \param c  input array 2 (complex)
   * \param d  input array 3 (real)
   */
   __host__ void realMulVConjVV(DeviceArray<cudaReal>& a, 
                                DeviceArray<cudaComplex> const & b,
                                DeviceArray<cudaComplex> const & c,
                                DeviceArray<cudaReal> const & d);
   
   /**
   * Performs qNew = (4 * (qr2 * expW2) - qr) / 3 elementwise, kernel wrapper
   * 
   * \param qNew  output array (a propagator slice)
   * \param qr  input array 1 (a propagator slice)
   * \param qr2  input array 2 (a propagator slice)
   * \param expW2  input array 3 (exp(-W[i]*ds/4) array)
   */
   __host__ void richardsonEx(DeviceArray<cudaReal>& qNew, 
                              DeviceArray<cudaReal> const & qr,
                              DeviceArray<cudaReal> const & qr2, 
                              DeviceArray<cudaReal> const & expW2);
   
   /**
   * Performs a[i] += b[i] * c[i] * d, kernel wrapper
   * 
   * \param a  output array
   * \param b  input array 1
   * \param c  input array 2
   * \param d  input scalar
   */
   __host__ void addEqMulVVc(DeviceArray<cudaReal>& a, 
                             DeviceArray<cudaReal> const & b,
                             DeviceArray<cudaReal> const & c, 
                             cudaReal const d);

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl<Propagator<D>>. A BlockTmpl<Propagator<D>>
   * has two Propagator<D> members and is derived from BlockDescriptor.
   *
   * \ingroup Rpg_Solvers_Module
   */
   template <int D>
   class Block : public BlockTmpl< Propagator<D> >
   {

   public:

      /**
      * Constructor.
      */
      Block();

      /**
      * Destructor.
      */
      ~Block();

      /**
      * Initialize discretization and allocate required memory.
      *
      * \param ds  desired (optimal) value for contour length step
      * \param mesh  Mesh<D> object - spatial discretization mesh
      * \param fft  FFT<D> object - Fourier transforms
      * \param useBatchedFFT  Flag indicating whether to use batched FFTs
      */
      void setDiscretization(double ds,
                             Mesh<D> const & mesh,
                             FFT<D> const & fft, 
                             bool useBatchedFFT = true);

      /**
      * Setup parameters that depend on the unit cell.
      *
      * \param unitCell  unit cell, defining cell dimensions (input)
      * \param waveList  container for properties of wavevectors (input)
      */
      void setupUnitCell(UnitCell<D> const & unitCell,
                         WaveList<D> const & waveList);

      /**
      * Setup parameters that depend on the unit cell.
      *
      * \param unitCell  unit cell, defining cell dimensions (input)
      */
      void setupUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set or reset block length.
      *
      * \param newLength  new block length
      */
      void setLength(double newLength);

      /**
      * Set or reset monomer statistical segment length.
      *
      * \param kuhn  new monomer statistical segment length.
      */
      void setKuhn(double kuhn);

      /**
      * Set solver for this block.
      *
      * \param w  chemical potential field for this monomer type
      */
      void setupSolver(RField<D> const & w);

      /**
      * Compute step of integration loop, from i to i+1.
      *
      * \param q  pointer to slice i of propagator q (input)
      * \param qNew pointer to slice i+1 of propagator q (output)
      */
      void step(RField<D> const & q, RField<D>& qNew);

      /**
      * Compute unnormalized concentration for block by integration.
      *
      * Upon return, grid point r of array cField() contains the
      * integral int ds q(r,s)q^{*}(r,L-s) times the prefactor,
      * where q(r,s) is the solution obtained from propagator(0),
      * and q^{*} is the solution of propagator(1),  and s is
      * a contour variable that is integrated over the domain
      * 0 < s < length(), where length() is the block length.
      *
      * \param prefactor  constant prefactor multiplying integral
      */
      void computeConcentration(double prefactor);

      /**
      * Compute derivatives of free energy w/ respect to cell parameters.
      *
      * The prefactor is the same as that used in computeConcentration.
      *
      * \param waveList  container for properties of wavevectors
      * \param prefactor  constant prefactor multiplying integral
      */
      void computeStress(WaveList<D> const & waveList, double prefactor);

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * \param n  unit cell parameter index
      */
      double stress(int n);

      /**
      * Return associated spatial Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Return associated FFT<D> object by const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Contour length step size.
      */
      double ds() const;

      /**
      * Number of contour length steps.
      */
      int ns() const;

      // Functions with non-dependent names from BlockTmpl<Propagator<D>>
      using BlockTmpl< Pscf::Rpg::Propagator<D> >::setKuhn;
      using BlockTmpl< Pscf::Rpg::Propagator<D> >::propagator;
      using BlockTmpl< Pscf::Rpg::Propagator<D> >::cField;
      using BlockTmpl< Pscf::Rpg::Propagator<D> >::length;
      using BlockTmpl< Pscf::Rpg::Propagator<D> >::kuhn;

      // Functions with non-dependent names from BlockDescriptor
      using BlockDescriptor::setId;
      using BlockDescriptor::setVertexIds;
      using BlockDescriptor::setMonomerId;
      using BlockDescriptor::setLength;
      using BlockDescriptor::id;
      using BlockDescriptor::monomerId;
      using BlockDescriptor::vertexIds;
      using BlockDescriptor::vertexId;
      using BlockDescriptor::length;

   private:

      /// Object to perform batched FFTs on the entire propagator at once
      FFTBatched<D> fftBatchedAll_;

      /// Object to perform batched FFTs on two fields at once
      FFTBatched<D> fftBatchedPair_;

      /// Stress contribution from this block
      FSArray<double, 6> stress_;

      /// exp(-K^2 b^2 ds/6) array on k-grid
      RField<D> expKsq_;

      /// exp(-K^2 b^2 ds/12) array on k-grid
      RField<D> expKsq2_;

      /// exp(-W[i] ds/2) array on r-grid
      RField<D> expW_;

      /// exp(-W[i] ds/4) array on r-grid
      RField<D> expW2_;

      /**
      * Workspace array containing two r-grid fields, stored on the device.
      * 
      * These workspace fields are stored contiguously in a single array to 
      * allow batched FFTs to be performed on both fields simultaneously, 
      * which occurs in step().
      */
      DeviceArray<cudaReal> qrPair_;

      /**
      * Workspace array containing two k-grid fields, stored on the device.
      * 
      * These workspace fields are stored contiguously in a single array to 
      * allow batched FFTs to be performed on both fields simultaneously, 
      * which occurs in step().
      */
      DeviceArray<cudaComplex> qkPair_;

      /// Container to store batched FFTs of q in contiguous memory
      DeviceArray<cudaComplex> qkBatched_;

      /// Container to store batched FFTs of q in contiguous memory
      DeviceArray<cudaComplex> qk2Batched_;

      // Propagators on r-grid
      RField<D> q1_;
      RField<D> q2_;

      /// Array of elements containing exp(-K^2 b^2 ds/6) on k-grid on host
      HostDArray<cudaReal> expKsq_h_;

      /// Array of elements containing exp(-K^2 b^2 ds/6) on k-grid on host
      HostDArray<cudaReal> expKsq2_h_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Pointer to associated FFT<D> object.
      FFT<D> const * fftPtr_;

      /// Pointer to associated UnitCell<D> object.
      UnitCell<D> const* unitCellPtr_;

      /// Pointer to associated WaveList<D> object.
      WaveList<D> const * waveListPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in discrete Fourier transform (DFT) k-grid
      int kSize_;

      /// Contour length step size (actual step size for this block)
      double ds_;

      /// Contour length step size (value input in param file)
      double dsTarget_;

      /// Number of chain contour positions (= # contour steps + 1).
      int ns_;

      /// Have arrays been allocated in setDiscretization ?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      /// Use batched FFTs to compute stress? (faster, but doubles memory use)
      bool useBatchedFFT_;

      /// Get associated UnitCell<D> as const reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Get the Wavelist as const reference
      WaveList<D> const & wavelist() const
      {  return *waveListPtr_; }

      /// Number of unit cell parameters
      int nParams_;

      /// Compute expKSq_ arrays.
      void computeExpKsq();

   };

   // Inline member functions

   // Get number of contour steps.
   template <int D>
   inline int Block<D>::ns() const
   {  return ns_; }

   // Get contour length step size.
   template <int D>
   inline double Block<D>::ds() const
   {  return ds_; }

   // Get derivative of free energy w/ respect to a unit cell parameter.
   template <int D>
   inline double Block<D>::stress(int n)
   {  return stress_[n]; }

   // Get Mesh by reference.
   template <int D>
   inline Mesh<D> const & Block<D>::mesh() const
   {
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   // Get FFT by reference.
   template <int D>
   inline FFT<D> const & Block<D>::fft() const
   {
      UTIL_ASSERT(fftPtr_);
      return *fftPtr_;
   }

   #ifndef RPG_BLOCK_TPP
   // Suppress implicit instantiation
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;
   #endif

}
}
#endif
