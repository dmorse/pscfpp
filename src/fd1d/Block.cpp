/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"

namespace Pscf { 
namespace Fd1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   Block::Block()
    : xMin_(0.0),
      xMax_(0.0),
      dx_(0.0),
      ds_(0.0),
      ns_(0),
      nx_(0)
   {
      propagator(0).setBlock(*this);
      propagator(1).setBlock(*this);
   }

   /*
   * Destructor.
   */
   Block::~Block()
   {}

   void Block::setGrid(double xMin, double xMax, int nx, int ns)
   {  
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      ns_ = ns;
      dx_ = (xMax - xMin)/double(nx_ - 1);
      ds_ = length()/double(ns_ - 1);
      dA_.allocate(nx_);
      uA_.allocate(nx_ - 1);
      dB_.allocate(nx_);
      uB_.allocate(nx_ - 1);
      v_.allocate(nx_);
      solver_.allocate(nx_);
      propagator(0).allocate(ns_, nx_);
      propagator(1).allocate(ns_, nx_);
      cField().allocate(nx_);
   }

   /*
   * Setup the contour length step algorithm.
   *
   * This implementation uses the Crank-Nicholson algorithm for stepping
   * the modified diffusion equation. One step of this algorithm, which
   * is implemented by the step() function, solves a matrix equation of 
   * the form
   *
   *         A q(i) = B q(i-1)
   *
   * where A and B are nx_ x nx_ symmetric tridiagonal matrices given by
   * 
   *           A = 1 + 0.5*ds_*H
   *           B = 1 + 0.5*ds_*H
   *
   * in which ds_ is the contour step and 
   *
   *           H = -(b^2/6)d^2/dx^2 + w 
   *
   * is a finite difference representation of the "Hamiltonian" operator,
   * in which b = kuhn() is the statistical segment length. 
   *
   * This function sets up arrays containing diagonal and off-diagonal 
   * elements of the matrices A and B, and computes the LU decomposition 
   * of matrix A. Arrays of nx_ diagonal elements of A and B are denoted 
   * by dA_ and dB_, respectively, while arrays of nx_ - 1 off-diagonal 
   * ("upper") elements of A and B are denoted by uA_ and uB_.
   */
   void Block::setupSolver(Block::WField const& w)
   {
      // Preconditions
      UTIL_CHECK(nx_ > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(dA_.capacity() == nx_);
      UTIL_CHECK(dB_.capacity() == nx_);
      UTIL_CHECK(uA_.capacity() == nx_ -1);
      UTIL_CHECK(uB_.capacity() == nx_ -1);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      // Initialize diagonals with identity and potential terms
      double halfdS = 0.5*ds_;
      for (int i = 0; i < nx_; ++i) {
         dA_[i] = 1.0 + halfdS*w[i];
         dB_[i] = 1.0 - halfdS*w[i];
      }

      // Add diagonal second derivative terms
      double bx = kuhn()/dx_;
      double c1 = bx*bx*ds_/12.0;
      double c2 = 2.0*c1;
      dA_[0] += c1;
      dB_[0] -= c1;
      for (int i = 1; i < nx_ - 1; ++i) {
         dA_[i] += c2;
         dB_[i] -= c2;
      }
      dA_[nx_ - 1] += c1;
      dB_[nx_ - 1] -= c1;

      // Assign off-diagonal second derivative terms
      // (Opposite sign of corresponding diagonals)
      for (int i = 0; i < nx_ - 1; ++i) {
         uA_[i] = -c1;
         uB_[i] = +c1;
      }

      // Compute the LU decomposition of the A matrix
      solver_.computeLU(dA_, uA_);

   }

   /*
   * Propagate solution by one step.
   *
   * This algorithm implements one step of the Crank-Nicholson algorithm.
   * To do so, it solves A q(i+1) = B q(i), where A and B are constant 
   * matrices defined in the documentation of the setupStep() function.
   */
   void Block::step(const QField& q, QField& qNew)
   {
      v_[0] = dB_[0]*q[0] + uB_[0]*q[1];
      for (int i = 1; i < nx_ - 1; ++i) {
         v_[i] = dB_[i]*q[i] + uB_[i-1]*q[i-1] + uB_[i]*q[i+1];
      }
      v_[nx_ - 1] = dB_[nx_-1]*q[nx_-1] + uB_[nx_-2]*q[nx_-2];
      solver_.solve(v_, qNew);
   }


   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Block::computeConcentration(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(nx_ > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(dA_.capacity() == nx_);
      UTIL_CHECK(dB_.capacity() == nx_);
      UTIL_CHECK(uA_.capacity() == nx_ -1);
      UTIL_CHECK(uB_.capacity() == nx_ -1);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == nx_) 

      // Initialize to zero at all points
      int i;
      for (i = 0; i < nx_; ++i) {
         cField()[i] = 0.0;
      }

      Propagator const & p0 = propagator(0);
      Propagator const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (i = 0; i < nx_; ++i) {
         cField()[i] += 0.5*p0.q(0)[i]*p1.q(ns_ - 1)[i];
      }
      for (int j = 1; j < ns_ - 1; ++j) {
         for (i = 0; i < nx_; ++i) {
            cField()[i] += p0.q(j)[i]*p1.q(ns_ - 1 - j)[i];
         }
      }
      for (i = 0; i < nx_; ++i) {
         cField()[i] += 0.5*p0.q(ns_ - 1)[i]*p1.q(0)[i];
      }

      // Normalize
      prefactor *= ds_;
      for (i = 0; i < nx_; ++i) {
         cField()[i] *= prefactor;
      }

   }

}
}
