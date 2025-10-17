#ifndef RPC_BLOCK_H
#define RPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <rpc/solvers/Propagator.h>       // template argument

#include <prdc/cpu/RField.h>              // member
#include <prdc/cpu/RFieldDft.h>           // member
#include <util/containers/DMatrix.h>      // member template
#include <util/containers/FSArray.h>      // member template

// Forward declarations
namespace Pscf {

   template <int D> class Mesh;
   namespace Prdc{
      template <int D> class UnitCell;
      namespace Cpu {
         template <int D> class FFT;
         template <int D> class WaveList;
      }
   }
   namespace Rpc{
      template <int D> class FieldIo;
   }

   // Explicit instantiation declarations for base classes
   extern template 
   class BlockTmpl< Rpc::Propagator<1>, Prdc::Cpu::RField<1> >;
   extern template 
   class BlockTmpl< Rpc::Propagator<2>, Prdc::Cpu::RField<2> >;
   extern template 
   class BlockTmpl< Rpc::Propagator<3>, Prdc::Cpu::RField<3> >;

}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Block within a linear or branched block polymer.
   *
   * A Block has two Propagator<D> members, and a RField<D> concentration 
   * field.
   *
   * \ref user_param_block_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Block : public BlockTmpl< Propagator<D>, RField<D> >
   {

   public:

      // Public type name aliases

      /// Base class.
      using Base = BlockTmpl< Propagator<D>, RField<D> >;

      /// Propagator type.
      using PropagatorT = Propagator<D>;

      /// Fast Fourier Transform type.
      using FFTT = FFT<D>;

      /// WaveList type.
      using WaveListT = WaveList<D>;

      /// FieldIo type.
      using FieldIoT = FieldIo<D>;

      // Public member functions

      /**
      * Constructor.
      */
      Block();

      /**
      * Destructor.
      */
      ~Block();

      /**
      * Create permanent associations with related objects.
      *
      * This function creates associations of this block with the Mesh, 
      * FFT, UnitCell and WaveList objects by storing their addresses.
      * It must be called before allocate().
      *
      * \param mesh  Mesh<D> object, spatial discretization meth
      * \param fft  FFT<D> object, Fast Fourier Transform 
      * \param cell  UnitCell<D> object, crystallographic unit cell
      * \param wavelist  WaveList<D>, container for wavevector properties
      */
      void associate(Mesh<D> const& mesh, 
                     FFT<D> const& fft, 
                     UnitCell<D> const& cell,
                     WaveList<D>& wavelist);

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
      void setupSolver(RField<D> const & w);

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
      void stepThread(RField<D> const & qin, RField<D>& qout) const;

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
      void stepBead(RField<D> const & qin, RField<D>& qout) const;

      /**
      * Apply the exponential field operator for the bead model. 
      *
      * This function applies exp( -w(r) ), where w(r) is the w-field for
      * the monomer type of this block. 
      *
      * \param q  slice of propagator q, modified in place
      */
      void stepFieldBead(RField<D> & q) const;

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
      void stepBondBead(RField<D> const & qin, RField<D>& qout) const;

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
      void stepHalfBondBead(RField<D> const & qin, RField<D>& qout) const;

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

      // Inherited public member functions with non-dependent names

      using Base::setKuhn;
      using Base::propagator;
      using Base::cField;
      using Base::kuhn;

      using Edge::setId;
      using Edge::setVertexIds;
      using Edge::setMonomerId;
      using Edge::setLength;
      using Edge::id;
      using Edge::monomerId;
      using Edge::vertexIds;
      using Edge::vertexId;
      using Edge::length;
      using Edge::nBead;

   private:

      // Private member data

      // In bead model, ds=1 by definition.

      /// Stress arising from this block.
      FSArray<double, 6> stress_;

      /// Array of elements containing exp(-K^2 b^2 ds/6)
      RField<D> expKsq_;

      /// Array containing exp(-W[i] ds/2) (thread) or exp(-W[i]) (bead)
      RField<D> expW_;

      /// Array of elements containing exp(-K^2 b^2 ds/(6*2)) 
      RField<D> expKsq2_;

      /// Array of elements containing exp(-W[i] (ds/2)*0.5) (thread model)
      RField<D> expW2_;

      /// Array of elements containing exp(+W[i]) (bead model)
      RField<D> expWInv_;

      /// Work array for real-space q field (step size ds)
      mutable RField<D> qr_;

      /// Work array for real-space q field (step size ds/2, thread model)
      mutable RField<D> qr2_;

      /// Work array for wavevector space field (step size ds)
      mutable RFieldDft<D> qk_;

      /// Work array for wavevector space field (step size ds/2)
      mutable RFieldDft<D> qk2_;

      /// Pointer to associated Mesh<D> object
      Mesh<D> const * meshPtr_;

      /// Pointer to associated FFT<D> object
      FFT<D> const * fftPtr_;

      /// Pointer to associated UnitCell<D> object
      UnitCell<D> const * unitCellPtr_;

      /// Pointer to associated WaveList<D> object (non-const)
      WaveList<D> * waveListPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh 
      int kSize_;

      /// Contour length step size (actual step size for this block)
      double ds_;

      /// Contour length step size (value input in param file)
      double dsTarget_;

      /// Number of contour grid points = # of contour steps + 1
      int ns_;

      /// Have arrays been allocated ?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      /// Get associated UnitCell<D> as const reference.
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /// Get associated WaveList<D> by const reference.
      WaveList<D> const & wavelist() const
      {  return *waveListPtr_; }

      /// Number of unit cell parameters.
      int nParams_;

      // Private member functions

      /**
      * Get associated spatial Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get associated FFT object by const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Compute expKSq arrays.
      */
      void computeExpKsq();

   };

   // Inline member functions

   // Get number of contour grid points, including end points.
   template <int D>
   inline int Block<D>::ns() const
   {  return ns_; }

   // Get contour step size.
   template <int D>
   inline double Block<D>::ds() const
   {  return ds_; }

   // Stress with respect to unit cell parameter n.
   template <int D>
   inline double Block<D>::stress(int n) const
   {  return stress_[n]; }

   // Get associated Mesh<D> object by const reference.
   template <int D>
   inline Mesh<D> const & Block<D>::mesh() const
   {
      UTIL_CHECK(meshPtr_);
      return *meshPtr_;
   }

   // Get associated FFT<D> object by const reference.
   template <int D>
   inline FFT<D> const & Block<D>::fft() const
   {
      UTIL_CHECK(fftPtr_);
      return * fftPtr_;
   }

   // Explicit instantiation declarations
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;

} // R1d
} // Pscf
#endif
