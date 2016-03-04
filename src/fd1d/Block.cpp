/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include "Domain.h"

namespace Pscf { 
namespace Fd1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   Block::Block()
    : domainPtr_(0),
      ds_(0.0),
      ns_(0)
   {
      propagator(0).setBlock(*this);
      propagator(1).setBlock(*this);
   }

   /*
   * Destructor.
   */
   Block::~Block()
   {}

   void Block::setDiscretization(Domain const & domain, double ds)
   {  
      UTIL_CHECK(length() > 0);

      // Set association to spatial domain
      domainPtr_ = &domain;

      // Set contour length discretization
      ns_ = floor(length()/ds + 0.5) + 1;
      if (ns_%2 == 0) {
         ns_ += 1;
      }
      ds_ = length()/double(ns_ - 1);

      // Allocate all required memory
      dA_.allocate(domain.nx());
      dB_.allocate(domain.nx());
      uA_.allocate(domain.nx() - 1);
      uB_.allocate(domain.nx() - 1);
      lA_.allocate(domain.nx() - 1);
      lB_.allocate(domain.nx() - 1);
      v_.allocate(domain.nx());
      solver_.allocate(domain.nx());
      propagator(0).allocate(ns_, domain.nx());
      propagator(1).allocate(ns_, domain.nx());
      cField().allocate(domain.nx());
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
   * where A and B are domain().nx() x domain().nx() symmetric tridiagonal 
   * matrices given by
   * 
   *           A = 1 + 0.5*ds_*H
   *           B = 1 + 0.5*ds_*H
   *
   * in which ds_ is the contour step and 
   *
   *           H = -(b^2/6)d^2/dx^2 + w 
   *
   * is a finite difference representation of the "Hamiltonian" 
   * operator, in which b = kuhn() is the statistical segment length. 
   *
   * This function sets up arrays containing diagonal and off-diagonal 
   * elements of the matrices A and B, and computes the LU 
   * decomposition of matrix A. Arrays of domain().nx() diagonal elements 
   * of A and B are denoted by dA_ and dB_, respectively, while arrays of 
   * domain().nx() - 1 upper and lower off-diagonal elements of A and B
   * are denoted by uA_, lA_, uB_, and lB_, respectively
   */
   void Block::setupSolver(Block::WField const& w)
   {
      // Preconditions
      UTIL_CHECK(domainPtr_);
      int nx = domain().nx();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(dA_.capacity() == nx);
      UTIL_CHECK(dB_.capacity() == nx);
      UTIL_CHECK(uA_.capacity() == nx - 1);
      UTIL_CHECK(uB_.capacity() == nx - 1);
      UTIL_CHECK(lA_.capacity() == nx - 1);
      UTIL_CHECK(lB_.capacity() == nx - 1);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      // Chemical potential terms in matrix A
      double halfdS = 0.5*ds_;
      for (int i = 0; i < nx; ++i) {
         dA_[i] = halfdS*w[i];
      }

      // Second derivative terms in matrix A
      double dx = domain().dx();
      double c0 = halfdS*kuhn()*kuhn()/(6.0*dx);
      double c1 = c0/dx;
      double c2 = 2.0*c1;
      dA_[0] += c2;
      uA_[0] = -c2;
      for (int i = 1; i < nx - 1; ++i) {
         dA_[i] += c2;
         uA_[i] = -c1;
         lA_[i-1] = -c1;
      }
      dA_[nx - 1] += c2;
      lA_[nx - 2] = -c2;

      // Additional first derivative terms (if any)
      GeometryMode mode = domain().geometryMode();
      if (mode != Planar) {
         double gamma = 0.0;
         if (mode == Cylindrical) {
            gamma = 1.0;
         } else 
         if (mode == Spherical) {
            gamma = 2.0;
         } else {
            UTIL_THROW("Unknown geometryMode");
         }
         double c3 = 0.5*gamma*c0;
         double xMin = domain().xMin();
         double x;
         for (int i = 1; i < nx - 1; ++i) {
            x = xMin + dx*i;
            uA_[i] -= c3/x;
            lA_[i-1] += c3/x;
         }
         if (xMin/dx < 0.1) {
            uA_[0] -= gamma*c2;
            dA_[0] += gamma*c2;
         }
      }

      // Construct matrix B
      for (int i = 0; i < nx; ++i) {
         dB_[i] = -dA_[i];
         dA_[i] += 1.0;
         dB_[i] += 1.0;
      }
      for (int i = 0; i < nx - 1; ++i) {
         uB_[i] = -uA_[i];
      }
      for (int i = 0; i < nx - 1; ++i) {
         lB_[i] = -lA_[i];
      }

      // Compute the LU decomposition of the A matrix
      solver_.computeLU(dA_, uA_, lA_);

   }

   /*
   * Propagate solution by one step.
   *
   * This function implements one step of the Crank-Nicholson algorithm.
   * To do so, it solves A q(i+1) = B q(i), where A and B are constant 
   * matrices defined in the documentation of the setupStep() function.
   */
   void Block::step(const QField& q, QField& qNew)
   {
      int nx = domain().nx();
      v_[0] = dB_[0]*q[0] + uB_[0]*q[1];
      for (int i = 1; i < nx - 1; ++i) {
         v_[i] = dB_[i]*q[i] + uB_[i-1]*q[i-1] + uB_[i]*q[i+1];
      }
      v_[nx - 1] = dB_[nx-1]*q[nx-1] + uB_[nx-2]*q[nx-2];
      solver_.solve(v_, qNew);
   }


   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Block::computeConcentration(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(domain().nx() > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(dA_.capacity() == domain().nx());
      UTIL_CHECK(dB_.capacity() == domain().nx());
      UTIL_CHECK(uA_.capacity() == domain().nx() -1);
      UTIL_CHECK(uB_.capacity() == domain().nx() -1);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == domain().nx()) 

      // Initialize cField to zero at all points
      int i;
      int nx = domain().nx();
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }

      Propagator const & p0 = propagator(0);
      Propagator const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(0)[i]*p1.q(ns_ - 1)[i];
      }
      for (int j = 1; j < ns_ - 1; ++j) {
         for (i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i]*p1.q(ns_ - 1 - j)[i];
         }
      }
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(ns_ - 1)[i]*p1.q(0)[i];
      }

      // Normalize
      prefactor *= ds_;
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }

   }

}
}
