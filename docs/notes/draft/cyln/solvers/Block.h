#ifndef CYLN_BLOCK_H
#define CYLN_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/BlockTmpl.h>       // base class template
#include "Propagator.h"                   // base class argument
#include <pscf/math/TridiagonalSolver.h>  // member
#include <cyln/field/Field.h>             // member
#include <cyln/field/FFT.h>               // member
#include <util/containers/DArray.h>       // member

namespace Pscf { 
namespace Cyln 
{ 

   class Domain;
   using namespace Util;

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl<Propagator>. A BlockTmpl<Propagator> has two 
   * Propagator members and is derived from BlockDescriptor.
   *
   * \ingroup Pscf_Cyln_Module
   */
   class Block : public BlockTmpl<Propagator>
   {

   public:

      /**
      * Monomer chemical potential field.
      */
      typedef Propagator::WField WField;

      /**
      * Constrained partition function q(r,s) for fixed s.
      */
      typedef Propagator::QField QField;

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
      * \param domain associated Domain object, with grid info
      * \param ds desired (optimal) value for contour length step
      */
      void setDiscretization(Domain const & domain, double ds);

      /**
      * Setup MDE solver for this block.
      */
      void setupSolver(WField const & w);

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
      * \param prefactor multiplying integral
      */ 
      void computeConcentration(double prefactor);

      /**
      * Compute step of integration loop, from i to i+1.
      */
      void step(QField const & q, QField& qNew);

      /**
      * Return associated domain by reference.
      */
      Domain const & domain() const;

      /**
      * Number of contour length steps.
      */
      int ns() const;

   private:

      // Data structures for pseudospectral algorithm
 
      // Fourier transform plan
      FFT fft_;

      // Work array for wavevector space field.
      Field<fftw_complex> qk_;

      // Array of elements containing exp(-W[i] ds/2)
      Field<double> expW_;

      // Array of elements containing exp(-K^2 b^2 ds/6)
      DArray<double> expKsq_;

      // Data structures for Crank-Nicholson algorithm
 
      /// Solver used in Crank-Nicholson algorithm
      TridiagonalSolver solver_;

      // Arrays dA_, uA_, lB_ dB_, uB_, luB_ contain elements of the 
      // the tridiagonal matrices A and B used in propagation by the
      // radial part of the Laplacian from step i to i + 1, which 
      // requires solution of a linear system A q(i+1) = B q(i).

      /// Diagonal elements of matrix A
      DArray<double> dA_;

      /// Off-diagonal upper elements of matrix A
      DArray<double> uA_;

      /// Off-diagonal lower elements of matrix A
      DArray<double> lA_;

      /// Diagonal elements of matrix B
      DArray<double> dB_;

      /// Off-diagonal upper elements of matrix B
      DArray<double> uB_;

      /// Off-diagonal lower elements of matrix B
      DArray<double> lB_;

      /// Work vector, dimension = nr = # elements in radial direction
      DArray<double> rWork_;

      /// Pointer to associated Domain object.
      Domain const * domainPtr_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /*
      * Setup data structures associated with Crank-Nicholson stepper
      * for radial part of Laplacian.
      */
      void setupRadialLaplacian(Domain& domain);

      /*
      * Setup data structures associated with pseudospectral stepper
      * for axial part of Laplacian, using FFTs.
      */
      void setupAxialLaplacian(Domain& domain);

   };

   // Inline member functions

   /// Get Domain by reference.
   inline Domain const & Block::domain() const
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   /// Get number of contour steps.
   inline int Block::ns() const
   {  return ns_; }

}
}
#endif
