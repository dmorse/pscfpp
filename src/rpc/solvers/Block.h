#ifndef RPC_BLOCK_H
#define RPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"                   // base class argument
#include <prdc/cpu/FFT.h>                 // member
#include <prdc/cpu/RField.h>              // member
#include <prdc/cpu/RFieldDft.h>           // member
#include <prdc/crystal/UnitCell.h>        // member
#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <pscf/mesh/Mesh.h>               // member
#include <util/containers/FArray.h>       // member template
#include <util/containers/DMatrix.h>      // member template

namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc{
     template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl< Propagator<D> >. A BlockTmpl< Propagator<D> >
   * has two Propagator<D> members and is derived from BlockDescriptor.
   *
   * \ref user_param_block_sec "Parameter File Format"
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
      * Initialize discretization and allocate required memory.
      *
      * This function creates associations of this block with the mesh, 
      * fft, and unit cell objects.
      *
      * \param mesh  spatial discretization mesh
      * \param fft  Fast Fourier Transform object
      * \param cell  unit cell object
      */
      void associate(Mesh<D> const& mesh, 
                     FFT<D> const& fft, 
                     UnitCell<D> const& cell);

      /**
      * Allocate memory and set contour step size.
      *
      * This function choses values for the number ns of contour 
      * variable grid points for this block and the associated step 
      * size length/(ns-1) for this block, and allocates memory for 
      * a variety of private arrays. 
      * 
      * The value for the number ns of contour variable grid points for 
      * this block so as to yield a value for the the actual step size 
      * length/(ns-1) as close as possible to the input parameter ds (the 
      * desired step size) consistent with the requirements that ns be an
      * odd integer and ns > 1. These requirements allow use of Simpson's 
      * rule for integration with respect to the contour variable s to
      * compute monomer concentration fields and stress contributions.
      *
      * \param ds desired (optimal) value for contour length step
      */
      void allocate(double ds);

      /**
      * Clear all internal data that depends on the unit cell parameters
      *
      * This function should be called once after every change in unit cell
      * parameters. The function marks all variables that depend on the 
      * unit cell parameters as being outdated and invalid. Such variables
      * are recomputed by lazy evaluation, just before they are needed.
      */
      void clearUnitCellData();

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
      * This should be called once after every change in w fields, the
      * unit cell parameters, block length or kuhn length, before
      * entering the loop used to solve the MDE for either propagator.
      * This function is called by Polymer<D>::compute.
      *
      * \param w chemical potential field for this monomer type
      */
      void setupSolver(RField<D> const & w);

      /**
      * Compute one step of solution of MDE for the thread model.
      *
      * This function is called internally by the PropagatorTmpl solve
      * function within a loop over steps. It is implemented in the
      * Block class because the same private data structures are needed
      * for the two propagators associated with a Block.
      *
      * \param qin  input slice of q, from step i
      * \param qout  output slice of q, from step i+1
      */
      void stepThread(RField<D> const & qin, RField<D>& qout);

      /**
      * Compute one step of solution of MDE for the bead model.
      *
      * This function is called internally by the PropagatorTmpl solve
      * function within a loop over steps. It is implemented in the
      * Block class because the same private data structures are needed
      * for the two propagators associated with a Block.
      *
      * \param qin  input slice of q, from step i
      * \param qout  output slice of q, for step i+1
      */
      void stepBead(RField<D> const & qin, RField<D>& qout);

      /**
      * Apply bond operator for the bead model. 
      *
      * \param qin  input slice of q, from step i
      * \param qout  ouptut slice of q, for step i+1
      */
      void stepBondBead(RField<D> const & qin, RField<D>& qout);

      /**
      * Apply the exponential field operator for the bead model. 
      *
      * \param q  slice of propagator q, modified in place.
      */
      void stepFieldBead(RField<D> & q);

      #if 0
      /**
      * Compute concentration (volume fraction) for block by integration.
      *
      * This should be called after both associated propagators are known.
      * This function is called by Polymer<D>::compute().
      *
      * Thread Model: 
      *
      * Bead Model:
      *
      *
      * \param prefactor  constant multiplying integral or sum
      */
      void computeConcentration(double prefactor);
      #endif

      /**
      * Compute the concentration for this block, for the thread model.
      *
      * This function is called by Polymer::compute if a thread model is
      * is used.
      *
      * The "prefactor" parameter must equal phi/(L q), where L is the 
      * total length of all blocks in the polymer species and q is the 
      * species partition function as
      *
      * Upon return, grid point r of array cField() contains the integal,
      * int ds q(r,s)q^{*}(r,L-s) times the prefactor parameter, where
      * q(r,s) and q^{*}(r,s) are propagators associated with different
      * directions, and the integral is taken over the length of the 
      * block. Simpson's rule is used for the integral.
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

      /**
      * Compute the spatial average of the product used to compute Q.
      *
      * This funciton computes the spatial average of the product 
      * q0[i]*q1[i], where q0 and q1 and are complementary propagator 
      * slices, and i is grid rank.
      */
      double averageProduct(RField<D> const& q0, RField<D> const& q1);

      /**
      * Compute the spatial average of the product used to compute Q.
      *
      * This computes the spatial average of the product 
      * q0[i]*q1[i]/exp(W[i]*ds), where q0 and q1 and are complementary 
      * propagator slices for a bead model, and i is grid rank. This is
      * used in the bead model for computation of Q from propagator 
      * slices associated with a bead that is owned by the propagator.
      */
      double averageProductBead(RField<D> const& q0, RField<D> const& q1);

      /**
      * Compute stress contribution for this block.
      *
      * This function is called by Polymer<D>::computeStress. The parameter
      * prefactor should be the same as that passed to the member function
      * computeConcentration.
      *
      * \param prefactor  constant multiplying integral over s
      */
      void computeStress(double prefactor);

      /**
      * Get associated spatial Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get associated FFT object by const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Get contour length step size.
      */
      double ds() const;

      /**
      * Get the number of contour length steps in this block.
      */
      int ns() const;

      /**
      * Get derivative of free energy w/ respect to unit cell parameter n.
      *
      * This function returns a value computed by a previous call to the
      * computeStress function.
      *
      * \param n index of unit cell parameter
      */
      double stress(int n) const;

      // Functions with non-dependent names from BlockTmpl< Propagator<D> >
      using BlockTmpl< Propagator<D> >::setKuhn;
      using BlockTmpl< Propagator<D> >::propagator;
      using BlockTmpl< Propagator<D> >::cField;
      using BlockTmpl< Propagator<D> >::length;
      using BlockTmpl< Propagator<D> >::kuhn;

      // Functions with non-dependent names from BlockDescriptor
      using BlockDescriptor::setId;
      using BlockDescriptor::setVertexIds;
      using BlockDescriptor::setMonomerId;
      using BlockDescriptor::setLength;
      using BlockDescriptor::id;
      using BlockDescriptor::monomerId;
      using BlockDescriptor::vertexIds;
      using BlockDescriptor::vertexId;
      using BlockDescriptor::ownsVertex;
      using BlockDescriptor::length;
      using BlockDescriptor::nBead;

   private:

      // Matrix to store derivatives of plane waves
      DMatrix<double> dGsq_;

      // Stress arising from this block
      FSArray<double, 6> stress_;

      // Fourier transform plan
      // FFT<D> fft_;

      // Array of elements containing exp(-K^2 b^2 ds/6)
      RField<D> expKsq_;

      // Array of elements containing exp(-W[i] ds/2) in thread model
      // or exp(-W[i] ds) in the bead model
      RField<D> expW_;

      // Array of elements containing exp(-K^2 b^2 ds/(6*2)) (thread model)
      RField<D> expKsq2_;

      // Array of elements containing exp(-W[i] (ds/2)*0.5) (thread model)
      RField<D> expW2_;

      // Array of elements containing exp(+W[i] ds) in bead model
      RField<D> expWInv_;

      // Work array for real-space field (step size ds)
      RField<D> qr_;

      // Work array for real-space field (step size ds/2, thread model)
      RField<D> qr2_;

      // Work array for wavevector space field (step size ds)
      RFieldDft<D> qk_;

      // Work array for wavevector space field (step size ds/2, thread model)
      RFieldDft<D> qk2_;

      // Pointer to associated Mesh<D> object
      Mesh<D> const* meshPtr_;

      // Pointer to associated FFT<D> object
      FFT<D> const* fftPtr_;

      // Pointer to associated UnitCell<D> object
      UnitCell<D> const* unitCellPtr_;

      // Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      // Number of wavevectors in wavevector mesh 
      int kSize_;

      // Contour length step size (actual step size for this block)
      double ds_;

      // Contour length step size (value input in param file)
      double dsTarget_;

      // Number of contour grid points = # of contour steps + 1
      int ns_;

      // Have arrays been allocated ?
      bool isAllocated_;

      // Are expKsq_ arrays up to date ? (initialize false)
      bool hasExpKsq_;

      /**
      * Access associated UnitCell<D> as reference.
      */
      UnitCell<D> const & unitCell() const
      {  return *unitCellPtr_; }

      /**
      * Compute dGsq_ matrix.
      */
      void computedGsq();

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

   /// Stress with respect to unit cell parameter n.
   template <int D>
   inline double Block<D>::stress(int n) const
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
      // return fft_;
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
