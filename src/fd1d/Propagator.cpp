/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"

namespace Fd1d
{

   using namespace Util;
   using namespace Chem;

   /*
   * Constructor.
   */
   Propagator::Propagator()
    : xMin_(0.0),
      xMax_(0.0),
      dx_(0.0),
      ds_(0.0),
      ns_(0),
      nx_(0)
   {}

   /*
   * Destructor.
   */
   Propagator::~Propagator()
   {}

   void Propagator::setGrid(double xMin, double xMax, int nx, int ns)
   {  
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      ns_ = ns;
      dx_ = (xMax - xMin)/double(nx_ - 1);
      ds_ = block().length()/double(ns_ - 1);
      qFields_.allocate(ns_);
      for (int i = 0; i < ns_; ++i) {
         qFields_[i].allocate(nx_);
      }
      dA_.allocate(nx_);
      uA_.allocate(nx_ - 1);
      dB_.allocate(nx_);
      uB_.allocate(nx_ - 1);
      v_.allocate(nx_);
      solver_.allocate(nx_);
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
   void Propagator::setupSolver(Propagator::WField const& w)
   {
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
   * Compute initial head QField from final tail QFields of sources.
   */
   void Propagator::computeHead()
   {

      // Reference to head of this propagator
      QField& qh = qFields_[0];

      // Initialize qh field to 1.0 at all grid points
      int ix;
      for (ix = 0; ix < nx_; ++ix) {
         qh[ix] = 1.0;
      }

      // Pointwise multiply tail QFields of all sources
      for (int is = 0; is < nSource(); ++is) {
         if (!source(is).isSolved()) {
            UTIL_THROW("Source not solved in computeHead");
         }
         QField const& qt = source(is).tail();
         for (ix = 0; ix < nx_; ++ix) {
            qh[ix] *= qt[ix];
         }
      }
   }

   /*
   * Propagate solution from step iStep to iStep + 1.
   *
   * This algorithm implements one step of the Crank-Nicholson algorithm.
   * To do so, it solves A q(i+1) = B q(i), where A and B are constant 
   * matrices defined in the documentation of the setupStep() function.
   */
   void Propagator::step(int iStep)
   {
      QField& q = qFields_[iStep];
      v_[0] = dB_[0]*q[0] + uB_[0]*q[1];
      for (int i = 1; i < nx_ - 1; ++i) {
         v_[i] = dB_[i]*q[i] + uB_[i-1]*q[i-1] + uB_[i]*q[i+1];
      }
      v_[nx_ - 1] = dB_[nx_-1]*q[nx_-1] + uB_[nx_-2]*q[nx_-2];
      solver_.solve(v_, qFields_[iStep + 1]);
   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   void Propagator::solve(const Propagator::WField& w)
   {
      computeHead();
      setupSolver(w);
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         step(iStep);
      }
      setIsSolved(true);
   }

   /*
   * Solve the modified diffusion equation with specified initial field.
   */
   void Propagator::solve(const Propagator::WField& w, const Propagator::QField& head) 
   {
      // Initialize head QField
      QField& qh = qFields_[0];
      for (int i = 0; i < nx_; ++i) {
         qh[i] = head[i];
      }

      // Setup solver and solve
      setupSolver(w);
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         step(iStep);
      }
      setIsSolved(true);
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Propagator::integrate(double prefactor, Propagator::CField& integral)
   {
      // Preconditions
      if (!hasPartner()) {
         UTIL_THROW("No partner");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner not solved");
      }
      if (integral.capacity() != nx_) {
         UTIL_THROW("integral not allocated or wrong size");
      }

      // Initialize to zero at all points
      int i;
      for (int i = 0; i < nx_; ++i) {
         integral[i] = 0.0;
      }

      // Evaluate unnormalized integral
      int j;
      for (j = 0; j < ns_; ++j) {
         for (i = 0; i < nx_; ++i) {
            integral[i] += qFields_[j][i]*partner().qFields_[ns_ - 1 - j][i];
         }
      }

      // Normalize
      prefactor *= dx_;
      for (i = 0; i < nx_; ++i) {
         integral[i] *= prefactor;
      }

   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   double Propagator::computeQ()
   {
      if (!isSolved()) {
         UTIL_THROW("Propagator is not solved.");
      }
      if (!hasPartner()) {
         UTIL_THROW("Propagator has no partner set.");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner propagator is not solved");
      }
      QField const& qh = head();
      QField const& qt = partner().tail();
      double Q = 0.0;
      for (int i = 0; i < nx_; ++i) {
         Q += qh[i]*qt[i];
      }
      return Q/double(nx_);
   }

}
