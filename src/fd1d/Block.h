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

      typedef typename Propagator::WField WField;
      typedef typename Propagator::QField QField;

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
      * \param xMin minimum value of space coordinate
      * \param xMax minimum value of space coordinate
      * \param nx number of spatial grid points
      * \param ns number of contour steps, # grid points - 1
      */
      void setGrid(double xMin, double xMax, int nx, int ns);

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
      * One step of integration loop, from i to i+1.
      */
      void step(QField const & q, QField& qNew);
 
   private:
 
      // Elements of tridiagonal matrices used in propagation
      DArray<double> dA_;
      DArray<double> dB_;
      DArray<double> uA_;
      DArray<double> uB_;

      // Work vector
      DArray<double> v_;

      // Solver used in Crank-Nicholson algorithm
      TridiagonalSolver solver_;

      /// Monomer statistical segment length.
      double step_;

      /// Minimum value of spatial coordinate.
      double xMin_;

      /// Maximum value of spatial coordinate.
      double xMax_;

      /// Spatial grid step size.
      double dx_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Number of spatial grid points.
      int nx_;

   };

} 
}
#endif
