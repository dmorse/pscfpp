#ifndef RPC_BLOCK_H
#define RPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <rpc/solvers/Propagator.h>       // base class argument

#include <prdc/cpu/RField.h>              // member
#include <prdc/cpu/RFieldDft.h>           // member
#include <util/containers/DMatrix.h>      // member template
#include <util/containers/FSArray.h>      // member template

// Forward references
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc{
     template <int D> class UnitCell;
     namespace Cpu {
        template <int D> class WaveList;
        template <int D> class FFT;
     }
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Block within a linear or branched block polymer.
   *
   * Derived from BlockTmpl< Propagator<D> >. A BlockTmpl< Propagator<D> >
   * has two Propagator<D> members, and is derived from class Pscf::Edge.
   *
   * \ref user_param_block_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
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
      * Create permanent associations with related objects.
      *
      * This function creates associations of this block with the mesh, 
      * fft, unit cell and wavelist objects by storing the addresses of 
      * these objects. This function must be called before the
      * allocate function.
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
      * For the bead model, if PolymerModel::isThread() is true, the
      * value of ns is given by nBead (the number of beads owned by the
      * block) plus one for each terminal vertex bead that this block 
      * does not own. 
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
      void stepThread(RField<D> const & qin, RField<D>& qout);

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
      void stepBead(RField<D> const & qin, RField<D>& qout);

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
      void stepBondBead(RField<D> const & qin, RField<D>& qout);

      /**
      * Apply the exponential field operator for the bead model. 
      *
      * This function applies exp( -w(r) ), where w(r) is the w-field for
      * the monomer type of this block. 
      *
      * \param q  slice of propagator q, modified in place
      */
      void stepFieldBead(RField<D> & q);

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
      * parameter. The sum is taken over all beads that are owned by this 
      * block. 
      *
      * \param prefactor  constant multiplying sum over beads
      */
      void computeConcentrationBead(double prefactor);

      /**
      * Compute the spatial average of a product used to compute Q.
      *
      * This function computes the spatial average of the product 
      * q0[i]*q1[i], where q0 and q1 and are complementary propagator 
      * slices, and i is a spatial mesh rank.
      */
      double averageProduct(RField<D> const& q0, RField<D> const& q1);

      /**
      * Compute the spatial average of a product used by the bead model.
      *
      * This computes the spatial average of the product 
      * q0[i]*q1[i]*exp(W[i]), where q0 and q1 and are complementary 
      * propagator slices for a bead model, and i is mesh rank. This 
      * is used in the bead model to compute Q from propagator slices 
      * associated with a bead that is owned by the propagator.
      */
      double averageProductBead(RField<D> const& q0, RField<D> const& q1);

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
      * Bead model: For the bead model, ns is equal to nBead plus the 
      * number of vertex beads that this block does not own, because the 
      * end-points of the mesh always correspond to the terminating beads, 
      * whichever block they are part of.
      * 
      * Thread mdoel: For the thread model, ns is always an odd number 
      * ns >= 3 chosen to give an even number of steps of length 
      * ds = length()/(ns() - 1) as close as possible to the target value 
      * of ds passed to the allocate function.
      *
      */
      int ns() const;

      // Functions with non-dependent names from BlockTmpl< Propagator<D> >
      using BlockTmpl< Propagator<D> >::setKuhn;
      using BlockTmpl< Propagator<D> >::propagator;
      using BlockTmpl< Propagator<D> >::cField;
      using BlockTmpl< Propagator<D> >::kuhn;

      // Functions with non-dependent names inherited from Edge
      using Edge::setId;
      using Edge::setVertexIds;
      using Edge::setMonomerId;
      using Edge::setLength;
      using Edge::id;
      using Edge::monomerId;
      using Edge::vertexIds;
      using Edge::vertexId;
      using Edge::ownsVertex;
      using Edge::length;
      using Edge::nBead;

   private:

      /// Stress arising from this block.
      FSArray<double, 6> stress_;

      /// Array of elements containing exp(-K^2 b^2 ds/6)
      RField<D> expKsq_;

      /// Array containing exp(-W[i] ds/2) (thread) or exp(-W[i] ds) (bead)
      RField<D> expW_;

      /// Array of elements containing exp(-K^2 b^2 ds/(6*2)) 
      RField<D> expKsq2_;

      /// Array of elements containing exp(-W[i] (ds/2)*0.5) (thread model)
      RField<D> expW2_;

      /// Array of elements containing exp(+W[i] ds) (bead model)
      RField<D> expWInv_;

      /// Work array for real-space field (step size ds)
      RField<D> qr_;

      /// Work array for real-space field (step size ds/2, thread model)
      RField<D> qr2_;

      /// Work array for wavevector space field (step size ds)
      RFieldDft<D> qk_;

      /// Work array for wavevector space field (step size ds/2)
      RFieldDft<D> qk2_;

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

   #ifndef RPC_BLOCK_TPP
   // Suppresse implicit instantiation
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;
   #endif

}
}
#endif
