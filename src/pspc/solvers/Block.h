#ifndef PSPC_BLOCK_H
#define PSPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"                   // base class argument
#include <pscf/solvers/BlockTmpl.h>       // base class template
#include <pscf/mesh/Mesh.h>               // member
#include <pscf/crystal/UnitCell.h>        // member
#include <pspc/field/RField.h>            // member
#include <pspc/field/RFieldDft.h>         // member
#include <pspc/field/FFT.h>               // member
#include <util/containers/FArray.h>       // member template
#include <util/containers/DMatrix.h>      // member template

namespace Pscf { 
   template <int D> class Mesh; 
   template <int D> class UnitCell;
}

namespace Pscf { 
namespace Pspc { 

   using namespace Util;

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl< Propagator<D> >. A BlockTmpl< Propagator<D> > 
   * has two Propagator<D> members and is derived from BlockDescriptor.
   *
   * \ref pscf_Block_page "Parameter File Format"
   *
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Block : public BlockTmpl< Propagator<D> >
   {

   public:

      // Typedefs

      /**
      * Generic field (r-grid format).
      */
      typedef typename Propagator<D>::Field Field;

      /**
      * Monomer chemical potential field (r-grid format).
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Constrained partition function q(r,s) for fixed s.
      */
      typedef typename Propagator<D>::QField QField;

      // Member functions

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
      * \param ds desired (optimal) value for contour length step
      * \param mesh spatial discretization mesh
      */
      void setDiscretization(double ds, const Mesh<D>& mesh);

      /**
      * Setup parameters that depend on the unit cell.
      *
      * This should be called once after every change in unit cell
      * parameters. Internal variables that depend on the unit cell
      * parameters are actually reset later, in in the setupSolver 
      * function just before solving the MDE, using a pointer to 
      * the unit cell that is set in this function. The associated
      * UnitCell<D> object should thus not be modified after calling
      * this function and before the MDEs are solved. 
      *
      * \param unitCell unit cell, defining cell dimensions
      */
      void setupUnitCell(const UnitCell<D>& unitCell);

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
      * This should be called once after every change in w fields, the
      * unit cell parameters, block length or kuhn lenght, before
      * entering the loop used to solve the MDE for either propagator.
      *
      * \param w chemical potential field for this monomer type
      */
      void setupSolver(WField const & w);

      /**
      * Compute one step of solution of MDE, from i to i+1.
      *
      * This function is called internally by the PropagatorTmpl solve
      * function within a loop over steps. It is implemented in the
      * Block class because the same private data structures are needed
      * for the two propagators associated with a Block.
      *
      * \param q  input value of QField, from step i
      * \param qNew  ouput value of QField, from step i+1
      */
      void step(QField const& q, QField& qNew);

      /**
      * Compute concentration (volume fraction) for block by integration.
      *
      * This should be called after both associated propagators are known.
      * Upon return, grid point r of array cField() contains the integal,
      * int ds q(r,s)q^{*}(r,L-s) times the prefactor parameter, where
      * q(r,s) is the solution obtained from propagator(0), q^{*}(r,s) is
      * the solution of propagator(1),  and s is a contour variable that 
      * is integrated over the domain 0 < s < length(), where length() 
      * is the block length. The "prefactor" parameter should be set to 
      * prefactor = phi/(L q), where phi is the overall volume fraction 
      * for this molecular species, L is the total number of monomers in 
      * the species, and q is the species partition function, i.e., the 
      * spatial average of q(r,L). This function is called by 
      * Polymer<D>::compute().
      *
      * \param prefactor  constant multiplying integral over s
      */ 
      void computeConcentration(double prefactor);

      /** 
      * Compute stress contribution for this block.
      *
      * This function is called by Polymer<D>::computeStress. The
      * prefactor parameter should be the same as that passed to 
      * function computeConcentration.   
      *   
      * \param prefactor  constant multiplying integral over s
      */  
      void computeStress(double prefactor);

      /**
      * Get associated spatial Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

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
      using BlockDescriptor::length;

   private:

      /// Matrix to store derivatives of plane waves 
      DMatrix<double> dGsq_;

      /// Stress arising from this block
      FSArray<double, 6> stress_;

      // Fourier transform plan
      FFT<D> fft_;

      // Array of elements containing exp(-K^2 b^2 ds/6)
      RField<D> expKsq_;

      // Array of elements containing exp(-W[i] ds/2)
      RField<D> expW_;

      // Array of elements containing exp(-K^2 b^2 ds/(6*2))
      RField<D> expKsq2_;

      // Array of elements containing exp(-W[i] (ds/2)*0.5)
      RField<D> expW2_;

      // Work array for real-space field.
      RField<D> qf_;

      // Work array for real-space field.
      RField<D> qr_;

      // Work array for real-space field.
      RField<D> qr2_;

      // Work array for wavevector space field.
      RFieldDft<D> qk_;

      // Work array for wavevector space field.
      RFieldDft<D> qk2_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const* meshPtr_;

      /// Pointer to associated UnitCell<D> object.
      UnitCell<D> const* unitCellPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform.
      IntVec<D> kMeshDimensions_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Have arrays been allocated in setDiscretization ?
      bool isAllocated_;

      /// Are expKsq_ arrays up to date ? (initialize false)
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

   #ifndef PSPC_BLOCK_TPP
   // Suppresse implicit instantiation
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;
   #endif

}
}
#endif
