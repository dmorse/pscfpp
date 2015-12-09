/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"

namespace Fd1d{

   using namespace Util;
   using namespace Chem;

   /*
   * Constructor.
   */
   Propagator::Propagator()
    : ds_(0.0),
      dx_(0.0),
      step_(0.0),
      ns_(0),
      nx_(0),
      hasHead_(false)
   {}

   void Propagator::init(int ns, int nx, double dx, double step)
   {
      ns_ = ns;
      nx_ = nx;
      dx_ = dx;
      step_ = step;
      ds_ = block().length()/double(ns_);
      qFields_.allocate(ns_ + 1);
      for (int i = 0; i <= ns_; ++i) {
         qFields_[i].allocate(nx_);
      }
      dA_.allocate(nx_);
      uA_.allocate(nx_ - 1);
      dB_.allocate(nx_);
      uB_.allocate(nx_ - 1);
      v_.allocate(nx_);
      solver_.allocate(nx_);
   }

   void Propagator::setHead(const Propagator::QField& head) 
   {
      QField& q = qFields_[0];
      for (int i = 0; i < nx_; ++i) {
         q[i] = head[i];
      }
      hasHead_ = true;
   }

   /*
   * Setup before the main propagation loop.
   *
   * The step routine propagates by one time step of the Crank-Nicholson 
   * algorithm, in which Aq(i) = Bq(i-1), with nx_ x nx_ matrices
   * 
   *           A = 1 + 0.5ds*H
   *           B = 1 + 0.5ds*H
   *           H = -(b^2/6)d^2/dx^2 + w 
   *
   * where H is the finite element representation of the "Hamiltonian" 
   * operator. This member function sets up the diagonal (dA_ and dB_) 
   * and off-diagonal (uA_ and uB_) elements of the tridiagonal matrices 
   * A and B, and computes and LU decomposition of A.
   */
   void Propagator::setup(Propagator::WField const& w)
   {
      // Initialize diagonals with identity and potential terms
      double halfdS = 0.5*ds_;
      for (int i = 0; i < nx_; ++i) {
         dA_[i] = 1.0 + halfdS*w[i];
         dB_[i] = 1.0 - halfdS*w[i];
      }

      // Add second derivative term
      double bx = step_/dx_;
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
      for (int i = 0; i < nx_ - 1; ++i) {
         uA_[i] -= c1;
         uB_[i] += c1;
      }
      solver_.computeLU(dA_, uA_);

      // Initialize head q-Field (if not already set externally)
      if (!hasHead_) {
         QField& q = qFields_[0];
         for (int j = 0; j < nx_; ++j) {
            q[j] = 1.0;
         }
         int j;
         for (int i = 0; i < nSource(); ++i) {
            QField const& qTail = source(i).tail();
            for (j = 0; j < nx_; ++j) {
               q[j] *= qTail[j];
            }
         }
      }

   }

   /*
   * Propagate from step iStep to iStep + 1.
   *
   * This algorithm solves A q(i+1) = B q(i), where A and B are matrices
   * defined in the documentation of the setup function.
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
      setup(w);
      for (int iStep = 0; iStep < ns_; ++iStep) {
         step(iStep);
      }
      hasHead_ = false;
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Propagator::integrate(Propagator::CField& integral)
   {
   }

}
