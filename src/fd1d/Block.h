#ifndef PSCF_BLOCK_H
#define PSCF_BLOCK_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include <pscf/BlockTmpl.h>
#include <pscf/TridiagonalSolver.h>    // member

namespace Pscf { 
namespace Fd1d 
{ 

   using namespace Util;

   class Block : public BlockTmpl<Propagator>
   {

   public:

      typedef Propagator::WField WField;
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
      * Initialize grid and allocate required memory.
      *
      * \param grid associated grid object
      * \param ds desired (optimal) value for contour length step
      */
      void setDiscretization(Grid const & grid, double ds);

      /**
      * Set Crank-Nicholson solver for this block.
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
 
   private:
 
      // Solver used in Crank-Nicholson algorithm
      TridiagonalSolver solver_;

      // Arrays dA_, uA_, dB_ and uB_ contain elements of the 
      // symmetric tridiagonal matrices A and B used in propagation 
      // from step i to i + 1, which requires solution of a linear 
      // system of the form: A q(i+1) = B q(i).

      // Diagonal elements of matrix A
      DArray<double> dA_;

      // Off-diagonal (upper or lower) elements of matrix A
      DArray<double> uA_;

      // Diagonal elements of matrix B
      DArray<double> dB_;

      // Off-diagonal elements of matrix B
      DArray<double> uB_;

      // Work vector
      DArray<double> v_;

      /// Pointer to associated Grid object.
      Grid const * gridPtr_;

      /// Monomer statistical segment length.
      double step_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Number of spatial grid points.
      int nx_;

      /// Return associated grid by reference.
      Grid const & grid() const;

   };

   // Inline member functions

   inline Grid const & Block::grid() const
   {   
      UTIL_ASSERT(gridPtr_);
      return *gridPtr_;
   }

} 
}
#endif
