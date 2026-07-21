#ifndef RPG_BLOCK_H
#define RPG_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/BlockTmpl.h>      // base class template

#include <pscf/cuda/CudaTp.h>            // template argument
#include <prdc/field/cuda/RField.h>            // member
#include <prdc/field/cuda/RFieldDft.h>         // member
#include <prdc/field/cuda/FFTBatched.h>        // member
#include <pscf/cuda/DeviceArray.h>       // member
#include <util/containers/FSArray.h>     // member

#include <prdc/field/cuda/resources.h>

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
      namespace Cuda {
         template <int D> class WaveList;
         template <int D> class FFT;
      }
   }
   namespace Rp {
      template <int D, class T> class Propagator;
      template <int D, class T> class Block;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Block within a branched polymer, for Cpp backend.
   *
   * A Block has two Propagator<D, CudaTp<D> > members and an 
   * RField<D, CudaTp<D> > monomer concentration field. 
   *
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Block<D, CudaTp<D> > 
     : public BlockTmpl< Rp::Propagator<D, CudaTp<D> >, RField<D, CudaTp<D> > >
   {

   public:

      // Public type name aliases

      /// Direct (parent) base class.
      using BlockTmplT = BlockTmpl< Rp::Propagator<D, CudaTp<D> >, RField<D, CudaTp<D> > > ;

      /// Propagator type (inherited).
      using typename BlockTmplT::PropagatorT;

      /// Field type (inherited).
      using typename BlockTmplT::FieldT;

      /// Fast Fourier Transform (FFT) type.
      using FFTT = FFT<D>;

      /// Wavelist type.
      using WaveListT = WaveList<D>;

      // Public member functions

      /**
      * Constructor.
      */
      Block();

      /**
      * Destructor.
      */
      ~Block();

      /// \name Setup (mutators)
      ///@{

      /**
      * Create permanent associations with related objects.
      *
      * This function creates associations of this block with the Mesh,
      * FFT, UnitCell, and WaveList objects by storing their addresses.
      * It must be called before allocate().
      *
      * \param mesh  Mesh<D> object - spatial discretization mesh
      * \param fft  FFT<D> object - Fourier transforms
      * \param cell  UnitCell<D> object - crystallographic unit cell
      * \param waveList  WaveList<D> object - properties of wavevectors
      */
      void associate(Mesh<D> const & mesh,
                     FFT<D> const & fft,
                     UnitCell<D> const & cell,
                     WaveList<D>& waveList);

      /**
      * Allocate memory and set contour step size for thread model.
      *
      * This function choses a value for the number ns of contour
      * variable grid points for this block, sets the step size, and
      * allocates memory for several private arrays. Spatial grid
      * dimensions are obtained from a pointers to the associated mesh.
      *
      * For the thread model, if PolymerModel::isThread() is true, the
      * value for the number ns of contour variable grid points for this
      * block is chosen to yield a value for the the actual step size
      * length/(ns-1) as close as possible to the input parameter ds (the
      * target step size), consistent with the requirements that ns be an
      * odd integer and ns > 1. These requirements allow use of Simpson's
      * rule for integration with respect to the contour variable s to
      * compute monomer concentration fields and stress contributions.
      *
      * For the bead model, if PolymerModel::isThread() is true, the value
      * of ns is given by nBead + 2.
      *
      * Precondition: associate() must be called before this function
      *
      * \param ds  target (optimal) value for contour length step size
      * \param useBatchedFFT  Flag indicating whether to use batched FFTs
      */
      void allocate(double ds, bool useBatchedFFT = true);

      /**
      * Clear all internal data that depends on lattice parameters.
      *
      * This method changes the internal hasExpKsq_ flag to false, so
      * that the expKsq arrays will need to be recalculated before
      * a step function can be called.
      */
      void clearUnitCellData();

      /**
      * Set or reset block length (only used in thread model).
      *
      * Precondition: PolymerModel::isThread()
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
      * This should be called once after every change in w fields, the
      * unit cell parameters, the block length or kuhn length, before
      * entering the loop used to solve the MDE for either propagator.
      * This function is called by Polymer<D>::compute.
      *
      * \param w  chemical potential field for this monomer type
      */
      void setupSolver(RField<D, CudaTp<D> > const & w);

      ///@}
      /// \name MDE Step Functions
      ///@{

      /**
      * Compute step of integration loop, from i to i+1.
      *
      * This function is called internally by the Propagator::solve
      * function within a loop over steps. It is implemented in the
      * Block class because the same private data structures are needed
      * for the two propagators associated with a Block.
      *
      * \param qin  slice i of propagator q (input)
      * \param qout  slice i+1 of propagator q (output)
      */
      void stepThread(RField<D, CudaTp<D> > const & qin, RField<D, CudaTp<D> >& qout) const;

      /**
      * Compute one step of solution of MDE for the bead model.
      *
      * This function is called internally by the Propagator::solve
      * function within a loop over steps. It is implemented in the
      * Block class because the same private data structures are needed
      * for the two propagators associated with a Block.
      *
      * \param qin  slice i of propagator q (input)
      * \param qout  slice i+1 of propagator q (output)
      */
      void stepBead(RField<D, CudaTp<D> > const & qin, RField<D, CudaTp<D> >& qout) const;

      /**
      * Apply the exponential field operator for the bead model.
      *
      * This function applies exp( -w(r)) in-place, where w(r) is the w-field
      * for the monomer type of this block.
      *
      * \param q  slice of propagator q, modified in placde
      */
      void stepFieldBead(RField<D, CudaTp<D> > & q) const;

      /**
      * Compute a bond operator for the bead model.
      *
      * This function applies exp( -|G|^2 b^2 / 6) in Fourier space.
      * To do so, it peforms a forward FFT, multplication in k-space,
      * and an inverse FFT.
      *
      * \param qin  slice i of propagator q (input)
      * \param qout  slice i+1 of propagator q (output)
      */
      void stepBondBead(RField<D, CudaTp<D> > const & qin, RField<D, CudaTp<D> >& qout) const;

      /**
      * Compute a half-bond operator for the bead model.
      *
      * This function applies exp( -|G|^2 b^2 / 12) in Fourier space.
      * To do so, it peforms a forward FFT, multplication in k-space,
      * and an inverse FFT.
      *
      * \param qin  slice i of propagator q (input)
      * \param qout  slice i+1 of propagator q (output)
      */
      void stepHalfBondBead(RField<D, CudaTp<D> > const & qin, RField<D, CudaTp<D> >& qout) const;

      ///@}
      /// \name Monomer Concentration Computation
      ///@{

      /**
      * Compute unnormalized concentration for block by integration.
      *
      * The "prefactor" parameter must equal phi/(L q), where L is the
      *  total length of all blocks in the polymer species and q is the
      * species partition function.
      *
      * Upon return, grid point r of array cField() contains the
      * integral int ds q(r,s)q^{*}(r,L-s) times the prefactor,
      * where q(r,s) is the solution obtained from propagator(0),
      * and q^{*} is the solution of propagator(1),  and s is
      * a contour variable that is integrated over the domain
      * 0 < s < length(), where length() is the block length.
      *
      * \param prefactor  constant multiplying integral over s
      */
      void computeConcentrationThread(double prefactor);

      /**
      * Compute the concentration for this block, using the bead model.
      *
      * This function is called by Polymer::compute if a bead model is
      * is used.
      *
      * The "prefactor" parameter must equal phi/(N q), where N is the
      * total number of beads owned by all blocks of the polymer, and
      * q is the species partition function.
      *
      * Upon return, grid point r of array cField() contains the sum
      * sum_s ds q(r,s) q^{*}(r,N-s) exp(W(r)*ds) over beads owned by
      * this block times the "prefactor" parameter, where q(r,s) and
      * q^{*}(r,s) are propagators associated with different directions.
      *
      * \param prefactor  constant multiplying sum over beads
      */
      void computeConcentrationBead(double prefactor);

      ///@}
      /// \name Stress Computation
      ///@{

      /**
      * Compute stress contribution for this block, using thread model.
      *
      * This function is called by Polymer<D>::computeStress. The
      * prefactor is equal to that used in computeConcentrationThread.
      *
      * \param prefactor  constant multiplying integral over s
      */
      void computeStressThread(double prefactor);

      /**
      * Compute stress contribution for this block, using bead model.
      *
      * This function is called by Polymer<D>::computeStress. The
      * prefactor is equal to that used in computeConcentrationBead.
      *
      * \param prefactor  constant multiplying sum over beads
      */
      void computeStressBead(double prefactor);

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * \param n  unit cell parameter index
      */
      double stress(int n) const;

      ///@}
      /// \name Accessors
      ///@{

      /**
      * Contour length step size.
      */
      double ds() const;

      /**
      * Number of contour length steps.
      */
      int ns() const;

      ///@}

      // Inherited functions with non-dependent names (selected)
      using BlockTmplT::propagator;
      using BlockTmplT::cField;

   private:

      /// Object to perform batched FFTs on the entire propagator at once
      FFTBatched<D> fftBatchedAll_;

      /// Object to perform batched FFTs on two fields at once
      FFTBatched<D> fftBatchedPair_;

      /// Stress contribution from this block
      FSArray<double, 6> stress_;

      /// Array containing exp(-K^2 b^2 ds / 6) for the thread model
      /// or exp(-K^2 b^2 /6) for the bead model
      RField<D, CudaTp<D> > expKsq_;

      /// Array containing exp(-W[i] ds/2) for the thread model
      /// or exp(-W[i]) for the bead model
      RField<D, CudaTp<D> > expW_;

      /// Array containing exp(-K^2 b^2 ds/12) (thread model)
      RField<D, CudaTp<D> > expKsq2_;

      /// Array containing exp(-W[i] ds/4) (thread model)
      RField<D, CudaTp<D> > expW2_;

      /// Array containing exp(+W[i]) (bead model)
      RField<D, CudaTp<D> > expWInv_;

      /**
      * Workspace array containing two r-grid fields, stored on the device.
      *
      * These workspace fields are stored contiguously in a single array to
      * allow batched FFTs to be performed on both fields simultaneously,
      * which occurs in stepThread().
      */
      mutable DeviceArray<cudaReal> qrPair_;

      /**
      * Workspace array containing two k-grid fields, stored on the device.
      *
      * These workspace fields are stored contiguously in a single array to
      * allow batched FFTs to be performed on both fields simultaneously,
      * which occurs in stepThread().
      */
      mutable DeviceArray<cudaComplex> qkPair_;

      // R-grid work space (used in productAverage)
      mutable RField<D, CudaTp<D> > qr_;

      // K-grid work space (used for FFT of q in stepBondBead)
      mutable RFieldDft<D> qk_;

      /// Container for batched FFTs of q0 (forward) in contiguous memory
      mutable DeviceArray<cudaComplex> q0kBatched_;

      /// Container for batched FFTs of q1 (reverse) in contiguous memory
      mutable DeviceArray<cudaComplex> q1kBatched_;

      // Slices of forward and reverse propagator on a k-grid (for stress)
      mutable RFieldDft<D> q0k_;
      mutable RFieldDft<D> q1k_;

      /// Const pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Const pointer to associated FFT<D> object.
      FFT<D> const * fftPtr_;

      /// Const pointer to associated UnitCell<D> object.
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated WaveList<D> object.
      WaveList<D> * waveListPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh
      int kSize_;

      /// Contour length step size (actual step size for this block)
      double ds_;

      /// Contour length step size (value input in param file)
      double dsTarget_;

      /// Number of chain contour positions (= # contour steps + 1).
      int ns_;

      /// Number of unit cell parameters
      int nParams_;

      /// Have arrays been allocated?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      /// Use batched FFTs to compute stress? (faster, but 2X memory use)
      bool useBatchedFFT_;

      /// Return associated spatial Mesh by const reference.
      Mesh<D> const & mesh() const;

      /// Return associated FFT<D> object by const reference.
      FFT<D> const & fft() const;

      /// Get associated UnitCell<D> by const reference (private).
      UnitCell<D> const & unitCell() const;

      /// Get the WaveList by non-const reference (private).
      WaveList<D> & waveList();

      /// Compute expKSq_ arrays.
      void computeExpKsq();

   };

   // Inline member functions

   // Get number of contour steps.
   template <int D> inline 
   int Block<D, CudaTp<D> >::ns() const
   {  return ns_; }

   // Get contour length step size.
   template <int D> inline 
   double Block<D, CudaTp<D> >::ds() const
   {  return ds_; }

   // Get derivative of free energy w/ respect to a unit cell parameter.
   template <int D> inline
   double Block<D, CudaTp<D> >::stress(int n) const
   {  return stress_[n]; }

   // Get Mesh by reference.
   template <int D> inline
   Mesh<D> const & Block<D, CudaTp<D> >::mesh() const
   {
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   // Get FFT by reference.
   template <int D> inline
   FFT<D> const & Block<D, CudaTp<D> >::fft() const
   {
      UTIL_ASSERT(fftPtr_);
      return *fftPtr_;
   }

   // Get associated UnitCell<D> by const reference (private).
   template <int D> inline
   UnitCell<D> const & Block<D, CudaTp<D> >::unitCell() const
   {  return *unitCellPtr_; }

   // Get the WaveList by non-const reference (private).
   template <int D> inline
   WaveList<D> & Block<D, CudaTp<D> >::waveList()
   {  return *waveListPtr_; }

} // namespace Rp
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {

   extern template
   class BlockTmpl< Rp::Propagator<1, CudaTp<1> >, Prdc::RField<1, CudaTp<1> > >;
   extern template
   class BlockTmpl< Rp::Propagator<2, CudaTp<2> >, Prdc::RField<2, CudaTp<2> > >;
   extern template
   class BlockTmpl< Rp::Propagator<3, CudaTp<3> >, Prdc::RField<3, CudaTp<3> > >;

   namespace Rp {
      extern template class Block<1, CudaTp<1> >;
      extern template class Block<2, CudaTp<2> >;
      extern template class Block<3, CudaTp<3> >;
   }
}
#endif
