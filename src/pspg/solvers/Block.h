#ifndef PSPG_BLOCK_H
#define PSPG_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"                   // base class argument
#include <prdc/cuda/RField.h>           // member
#include <prdc/cuda/RFieldDft.h>        // member
#include <prdc/cuda/FFT.h>               // member
#include <pspg/field/FFTBatched.h>        // member
#include <pspg/solvers/WaveList.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <util/containers/FArray.h>

namespace Pscf {

   template <int D> class Mesh;

namespace Pspg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl<Propagator<D>>. A BlockTmpl<Propagator<D>>
   * has two Propagator<D> members and is derived from BlockDescriptor.
   *
   * \ingroup Pspg_Solvers_Module
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
      */
      void setDiscretization(double ds,
                             Mesh<D> const & mesh,
                             FFT<D> const & fft);

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
      * \param length  new block length
      */
      void setLength(double length);

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
      * Initialize FFT and batch FFT classes.
      */
      void setupFFT();

      /**
      * Compute step of integration loop, from i to i+1.
      *
      * \param q  pointer to current slice of propagator q (input)
      * \param qNew pointer to current slice of propagator q (output)
      */
      void step(cudaReal const * q, cudaReal* qNew);

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
      using BlockTmpl< Pscf::Pspg::Propagator<D> >::setKuhn;
      using BlockTmpl< Pscf::Pspg::Propagator<D> >::propagator;
      using BlockTmpl< Pscf::Pspg::Propagator<D> >::cField;
      using BlockTmpl< Pscf::Pspg::Propagator<D> >::length;
      using BlockTmpl< Pscf::Pspg::Propagator<D> >::kuhn;

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

      /// Number of GPU blocks, set in setDiscretization
      int nBlocks_;

      /// Number of GPU threads per block, set in setDiscretization
      int nThreads_;

      /// Batched FFT, used in computeStress
      FFTBatched<D> fftBatched_;

      /// Stress conntribution from this block
      FArray<double, 6> stress_;

      // Array of elements containing exp(-K^2 b^2 ds/6) on k-grid
      RField<D> expKsq_;

      // Array of elements containing exp(-K^2 b^2 ds/12) on k-grid
      RField<D> expKsq2_;

      // Array of elements containing exp(-W[i] ds/2) on r-grid
      RField<D> expW_;

      // Array of elements containing exp(-W[i] ds/4) on r-grid
      RField<D> expW2_;

      // Work arrays for r-grid fields
      RField<D> qr_;
      RField<D> qr2_;

      // Work arrays for wavevector space (k-grid) field
      RFieldDft<D> qk_;
      RFieldDft<D> qk2_;

      // Batched FFTs of q
      cudaComplex* qkBatched_;
      cudaComplex* qk2Batched_;

      // Propagators on r-grid
      RField<D> q1_;
      RField<D> q2_;

      cudaReal* d_temp_;
      cudaReal* temp_;

      cudaReal* expKsq_host;
      cudaReal* expKsq2_host;

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

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # s nodes - 1.
      int ns_;

      /// Have arrays been allocated in setDiscretization ?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      /// Get associated UnitCell<D> as const reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Get the Wavelist as const reference
      WaveList<D> const & wavelist() const
      {  return *waveListPtr_; }

      /// Number of unit cell parameters
      int nParams_;

      /**
      * Compute expKSq_ arrays.
      */
      void computeExpKsq();

   };

   // Inline member functions

   /// Get number of contour steps.
   template <int D>
   inline int Block<D>::ns() const
   {  return ns_; }

   /// Get number of contour steps.
   template <int D>
   inline double Block<D>::ds() const
   {  return ds_; }

   /// Get derivative of free energy w/ respect to a unit cell parameter.
   template <int D>
   inline double Block<D>::stress(int n)
   {  return stress_[n]; }

   /// Get Mesh by reference.
   template <int D>
   inline Mesh<D> const & Block<D>::mesh() const
   {
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   /// Get FFT by reference.
   template <int D>
   inline FFT<D> const & Block<D>::fft() const
   {
      UTIL_ASSERT(fftPtr_);
      return *fftPtr_;
   }

   #ifndef PSPG_BLOCK_TPP
   // Suppress implicit instantiation
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;
   #endif

}
}
#endif
