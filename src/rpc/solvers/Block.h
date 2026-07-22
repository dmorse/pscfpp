#ifndef RPC_BLOCK_H
#define RPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <pscf/cpu/CppTp.h>           // template argument

#include <prdc/field/cpu/RField.h>        // member
#include <prdc/field/cpu/RFieldDft.h>     // member
#include <util/containers/FSArray.h>      // member

// Forward declarations
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc{
      template <int D> class UnitCell;
      template <int D, class T> class FFT;
      template <int D, class T> class WaveList;
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

   /**
   * Block within a linear or branched block polymer.
   *
   * A Block has two Propagator<D, Types<D> > members, and an 
   * RField<D, CppTp<D> > concentration field.
   *
   * \ref user_param_block_sec "Manual Page"
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Block<D, CppTp<D> > 
    : public BlockTmpl< Rp::Propagator<D, CppTp<D> >, RField<D, CppTp<D> > >
   {

   public:

      // Public type name aliases

      /// Direct (parent) base class.
      using BlockTmplT = BlockTmpl< Rp::Propagator<D, CppTp<D> >, RField<D, CppTp<D> > >;

      /// Propagator type (inherited).
      using typename BlockTmplT::PropagatorT;

      /// Field type (inherited).
      using typename BlockTmplT::FieldT;

      // Public member functions

      /**
      * Constructor.
      */
      Block();

      /**
      * Destructor.
      */
      ~Block();

      /// \name Initialization and State Mutators
      ///@{
     
      /**
      * Create permanent associations with related objects.
      *
      * This function creates associations of this block with the Mesh,
      * FFT, UnitCell and WaveList objects by storing their addresses.
      * It must be called before allocate().
      *
      * \param mesh  Mesh<D> object, spatial discretization meth
      * \param fft  FFT<D, CppTp<D> > object, Fast Fourier Transform
      * \param cell  UnitCell<D> object, crystallographic unit cell
      * \param waveList  WaveList<D, CppTp<D> >, container for wavevector properties
      */
      void associate(Mesh<D> const& mesh,
                     FFT<D, CppTp<D> > const& fft,
                     UnitCell<D> const& cell,
                     WaveList<D, CppTp<D> >& waveList);

      /**
      * Allocate memory and set contour step size.
      *
      * This function choses a value for the number ns of contour
      * variable grid points for this block, sets the step size, and
      * allocates memory for several private arrays. Spatial grid
      * dimensions are obtained from a pointers to the associated mesh.
      * The associate function must be called before this function.
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
      * The value of input parameter ds is ignored for the bead model.
      *
      * \param ds desired (optimal) value for contour length step
      */
      void allocate(double ds);

      /**
      * Clear all internal data that depends on the unit cell parameters
      *
      * This function must be called once after every time the unit cell
      * parameters change. The function marks all class member variables
      * that depend on the unit cell parameters as being outdated. All
      * such variables are then recomputed just before they are needed.
      */
      void clearUnitCellData();

      /**
      * Set or reset block length (only used in thread model).
      *
      * Precondition: PolymerModel::isThread(). An Exception is thrown if
      * this function is called when PolymerModel::isThread() is false.
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
      * Set up the MDE solver for this block.
      *
      * This should be called once after every change in w fields, before
      * entering the loop used to solve the MDE for either propagator.
      * This function is called by Polymer<D>::compute.
      *
      * \param w chemical potential field for this monomer type
      */
      void setupSolver(RField<D, CppTp<D> > const & w);

      ///@}
      /// \name MDE Step Functions
      ///@{
 
      /**
      * Compute one step of solution of MDE for the thread model.
      *
      * This function is called internally by the Propagator::solve
      * function within a loop over steps. It is implemented in the Block
      * class because the same private data structures are needed for the
      * two propagators associated with a Block.
      *
      * \param qin  input slice of q, from step i
      * \param qout  output slice of q, from step i+1
      */
      void stepThread(RField<D, CppTp<D> > const & qin, RField<D, CppTp<D> >& qout) const;

      /**
      * Compute one step of solution of MDE for the bead model.
      *
      * This function is called internally by the Propagator::solve
      * function within a loop over steps. It is implemented in the Block
      * class because the same private data structures are needed for the
      * two propagators associated with a Block.
      *
      * \param qin  input slice of q, from step i
      * \param qout  output slice of q, for step i+1
      */
      void stepBead(RField<D, CppTp<D> > const & qin, RField<D, CppTp<D> >& qout) const;

      /**
      * Apply the exponential field operator for the bead model.
      *
      * This function applies exp( -w(r) ), where w(r) is the w-field for
      * the monomer type of this block.
      *
      * \param q  slice of propagator q, modified in place
      */
      void stepFieldBead(RField<D, CppTp<D> > & q) const;

      /**
      * Apply a bond operator for the bead model.
      *
      * This function applies exp( nabla^2 b^2 / 6 ), where nabla^2
      * denotes a Laplacian operator with eigenvalues given by -G^2 for
      * reciprocal lattice vectors.
      *
      * \param qin  input slice of q, from step i
      * \param qout  ouptut slice of q, for step i+1
      */
      void stepBondBead(RField<D, CppTp<D> > const & qin, RField<D, CppTp<D> >& qout) const;

      /**
      * Apply a half-bond operator for the bead model.
      *
      * This function applies exp( nabla^2 b^2 / 12 ), where nabla^2
      * denotes a Laplacian operator with eigenvalues given by -G^2 for
      * reciprocal lattice vectors. It is used in the Propagator::solve
      * function to deal with half-bonds at block ends.
      *
      * \param qin  input slice of q, from step i
      * \param qout  ouptut slice of q, for step i+1
      */
      void stepHalfBondBead(RField<D, CppTp<D> > const & qin, RField<D, CppTp<D> >& qout) const;

      ///@}
      /// \name Monomer Concentration Computation
      ///@{
 
      /**
      * Compute the concentration for this block, for the thread model.
      *
      * This function is called by Polymer::compute if a thread model is
      * is used.
      *
      * The "prefactor" parameter must equal \f$ \phi / (L_{tot} Q) \f$,
      * where \f$ \phi \f$ is the species volume fraction, \f$ L_{tot} \f$
      * is the total length of all blocks in this polymer species and Q
      * is the species partition function.
      *
      * Upon return, grid point r of the array returned by the member
      * function cField() contains the integal
      * \f[
      *      p \int_{0}^{l} ds q_{0}(r,s) q_{1}(r, L - s)
      * \f]
      * where \f$ q_{0}(r,s) \f$ and \f$ q_{1}(r,s) \f$ are propagators
      * associated with different directions, \f$ p \f$ is the prefactor
      * parameter, and the integral is taken over the length \f$ L \f$
      * of this block. Simpson's rule is used for the integral with
      * respect to s.
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
      * The "prefactor" parameter must equal \f$ \phi /(N_{tot} Q) \f$,
      * where \f$ \phi \f$ is the species volume fraction, \f$ N_{tot} \f$
      * is the total number of beads in all blocks of the polymer, and
      * \f$ Q \f$ is the species partition function.
      *
      * Upon return, grid point r of the array returned by member function
      * cField() contains the sum
      * \f[
      *      p \sum_{s} q_{0}(r,s) q_{1}(r, N-s) \exp(W(r)*ds)
      * \f]
      * where \f$ q_{0}(r,s) \f$ and \f$ q_{1}(r, N-s) \f$ denote
      * complementary propagator slices associated with different
      * directions but the same bead, and \f$ p \f$ is the prefactor
      * parameter. The sum is taken over all beads in this block.
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
      * prefactor parameter must be equal to that passed to function
      * computeConcentrationThread(double ).
      *
      * \param prefactor  constant multiplying integral over s
      */
      void computeStressThread(double prefactor);

      /**
      * Compute stress contribution for this block, using bead model.
      *
      * This function is called by Polymer<D>::computeStress. The
      * prefactor parameter must be equal to that passed to function
      * computeConcentrationBead(double ).
      *
      * \param prefactor  constant multiplying sum over beads
      */
      void computeStressBead(double prefactor);

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * This function returns a value computed by a previous call to the
      * computeStress function.
      *
      * \param n  index of the unit cell parameter
      */
      double stress(int n) const;

      ///@}
      /// \name Accessors
      ///@{
  
      /**
      * Get contour length step size.
      */
      double ds() const;

      /**
      * Get the number of contour grid points, including end points.
      *
      * Thread mdoel: For the thread model, ns is always an odd number
      * ns >= 3 chosen to give an even number ns - 1 of steps of with
      * a step length ds = length/(ns - 1) as close as possible to the
      * target value of ds passed to the allocate function.
      *
      * Bead model: For the bead model, ns is equal to nBead + 2, so as
      * to include slices for all beads and two attached phantom vertices.
      * Beads are indexed 1, ..., ns - 2, vertices have indices 0 and
      * ns - 1.
      */
      int ns() const;

      // Inherited inherited functions with non-dependent names (selected)
      using BlockTmplT::propagator;
      using BlockTmplT::cField;

   private:

      // Private member data

      /// Stress arising from this block.
      FSArray<double, 6> stress_;

      /// Array of elements containing exp(-K^2 b^2 ds/6)
      RField<D, CppTp<D> > expKsq_;

      /// Array containing exp(-W[i] ds/2) (thread) or exp(-W[i]) (bead)
      RField<D, CppTp<D> > expW_;

      /// Array of elements containing exp(-K^2 b^2 ds/(6*2))
      RField<D, CppTp<D> > expKsq2_;

      /// Array of elements containing exp(-W[i] (ds/2)*0.5) (thread model)
      RField<D, CppTp<D> > expW2_;

      /// Array of elements containing exp(+W[i]) (bead model)
      RField<D, CppTp<D> > expWInv_;

      /// Work array for real-space q field (step size ds)
      mutable RField<D, CppTp<D> > qr_;

      /// Work array for real-space q field (step size ds/2, thread model)
      mutable RField<D, CppTp<D> > qr2_;

      /// Work array for wavevector space field (step size ds)
      mutable RFieldDft<D, CppTp<D> > qk_;

      /// Work array for wavevector space field (step size ds/2)
      mutable RFieldDft<D, CppTp<D> > qk2_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      /// Pointer to associated FFT<D, CppTp<D> > object
      FFT<D, CppTp<D> > const * fftPtr_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated WaveList<D, CppTp<D> > object (non-const)
      WaveList<D, CppTp<D> > * waveListPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh
      int kSize_;

      /// Contour step size (actual step size for this block)
      // In bead model, ds=1 by definition.
      double ds_;

      /// Target contour step size for thread model (from param file)
      double dsTarget_;

      /// Number of contour grid points = # of contour steps + 1
      int ns_;

      /// Number of unit cell parameters.
      int nParams_;

      /// Have arrays been allocated ?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      // Private member functions

      /// Get associated spatial Mesh by const reference (private).
      Mesh<D> const & mesh() const;

      /// Get associated FFT object by const reference (private).
      FFT<D, CppTp<D> > const & fft() const;

      /// Get associated UnitCell<D> as const reference (private).
      UnitCell<D> const & unitCell() const;

      /// Get associated WaveList<D, CppTp<D> > by non-const reference (private).
      WaveList<D, CppTp<D> > & waveList();

      /**
      * Compute expKSq arrays.
      */
      void computeExpKsq();

   };

   // Inline member functions

   // Get number of contour grid points, including end points.
   template <int D>
   inline int Block<D, CppTp<D> >::ns() const
   {  return ns_; }

   // Get contour step size.
   template <int D>
   inline double Block<D, CppTp<D> >::ds() const
   {  return ds_; }

   // Stress with respect to unit cell parameter n.
   template <int D>
   inline double Block<D, CppTp<D> >::stress(int n) const
   {  return stress_[n]; }

   // Private inline member function definitions

   // Get associated Mesh<D> object by const reference (private).
   template <int D>
   inline Mesh<D> const & Block<D, CppTp<D> >::mesh() const
   {
      UTIL_CHECK(meshPtr_);
      return *meshPtr_;
   }

   // Get associated FFT<D, CppTp<D> > object by const reference (private).
   template <int D>
   inline FFT<D, CppTp<D> > const & Block<D, CppTp<D> >::fft() const
   {
      UTIL_CHECK(fftPtr_);
      return * fftPtr_;
   }

   // Get associated UnitCell<D> by const reference (private).
   template <int D>
   UnitCell<D> const & Block<D, CppTp<D> >::unitCell() const
   {
      UTIL_CHECK(unitCellPtr_);
      return *unitCellPtr_;
   }

   // Get associated WaveList<D, CppTp<D> > by reference (private).
   template <int D>
   WaveList<D, CppTp<D> >& Block<D, CppTp<D> >::waveList()
   {
      UTIL_CHECK(waveListPtr_);
      return *waveListPtr_;
   }

} // namespace R1d
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {

   extern template
   class BlockTmpl< Rp::Propagator<1, CppTp<1> >, Prdc::RField<1, CppTp<1> > >;
   extern template
   class BlockTmpl< Rp::Propagator<2, CppTp<2> >, Prdc::RField<2, CppTp<2> > >;
   extern template
   class BlockTmpl< Rp::Propagator<3, CppTp<3> >, Prdc::RField<3, CppTp<3> > >;

   namespace Rp {
      extern template class Block<1, CppTp<1> >;
      extern template class Block<2, CppTp<2> >;
      extern template class Block<3, CppTp<3> >;
   }

}
#endif
