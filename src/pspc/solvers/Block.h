#ifndef PSPC_BLOCK_H
#define PSPC_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
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
#include <util/containers/DMatrix.h>

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
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Block : public BlockTmpl< Propagator<D> >
   {

   public:

      /**
      * Generic field (base class)
      */
      typedef typename Propagator<D>::Field Field;

      /**
      * Monomer chemical potential field.
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
      * \param unitCell unit cell, defining cell dimensions
      */
      void setupUnitCell(const UnitCell<D>& unitCell);

      /**
      * Set solver for this block.
      *
      * \param w  chemical potential field for this monomer type
      */
      void setupSolver(WField const & w);

      /**
      * Compute step of integration loop, from i to i+1.
      */
      void step(QField const & q, QField & qNew);

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
      * \param prefactor constant multiplying integral
      */ 
      void computeConcentration(double prefactor);

      /** 
      * Compute stress contribution for this block.
      *
      * The prefactor should be the same as that used in function
      * computeConcentration.   
      *   
      * \param prefactor constant multiplying integral
      */  
      void computeStress(double prefactor);

      /**
      * Return associated spatial Mesh by reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Contour length step size.
      */
      double ds() const;

      /**
      * Number of contour length steps.
      */
      int ns() const;

      /**
      * Get derivative of free energy w/ respect to unit cell parameter n.
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

      /// Pointer to associated Mesh<D>
      //const Mesh<D>* meshPtr_;

      /// Matrix to store derivatives of plane waves 
      DMatrix<double> dGsq_;

      /**
      * Compute dGsq_.
      */
      void computedGsq();

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

      /// Pointer to associated UnitCell<D>
      UnitCell<D> const* unitCellPtr_;

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /** 
      * Access associated UnitCell<D> as reference.
      */  
      UnitCell<D> const & unitCell() const { return *unitCellPtr_; }

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
   extern template class Block<1>;
   extern template class Block<2>;
   extern template class Block<3>;
   #endif
}
}
#endif
