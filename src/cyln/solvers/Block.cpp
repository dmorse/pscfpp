/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <cyln/misc/Domain.h>

namespace Pscf { 
namespace Cyln
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
      UTIL_CHECK(domain.nr() > 1);
      UTIL_CHECK(domain.nz() > 1);
      UTIL_CHECK(ds > 0.0);

      // Set association to spatial domain
      domainPtr_ = &domain;

      // Set contour length discretization
      ns_ = floor(length()/ds + 0.5) + 1;
      if (ns_%2 == 0) {
         ns_ += 1;
      }
      ds_ = length()/double(ns_ - 1);

      int nr = domain.nr();
      int nz = domain.nz();

      // Allocate propagators and cField
      propagator(0).allocate(ns_, nr, nz);
      propagator(1).allocate(ns_, nr, nz);
      cField().allocate(nr, nz);

      // Allocate work arrays for radial Crank-Nicholson
      dA_.allocate(nr);
      dB_.allocate(nr);
      uA_.allocate(nr - 1);
      uB_.allocate(nr - 1);
      lA_.allocate(nr - 1);
      lB_.allocate(nr - 1);
      v_.allocate(nr);
      solver_.allocate(nr);

   }

   /*
   * Setup the Crank-Nicholson work arrays.
   *
   * This implementation uses the Crank-Nicholson algorithm for stepping
   * the modified diffusion equation. One step of this algorithm, which
   * is implemented by the step() function, solves a matrix equation of 
   * the form
   *
   *         A q(i) = B q(i-1)
   *
   * where A and B are nr x nr symmetric tridiagonal matrices given by
   * 
   *           A = 1 + 0.5*ds_*H
   *           B = 1 - 0.5*ds_*H
   *
   * in which ds_ is the contour step and 
   *
   *           H = -(b^2/6)d^2/dr^2 
   *
   * is a finite difference representation of the radial part of the
   * "Hamiltonian" operator, in which b = kuhn() is the statistical 
   * segment length. 
   *
   * This function sets up arrays containing diagonal and off-diagonal 
   * elements of the matrices A and B, and computes the LU 
   * decomposition of matrix A. Arrays of the nr diagonal elements of
   * A and B are denoted by dA_ and dB_, respectively, while arrays 
   * of nr - 1 upper and lower off-diagonal elements of A and B are
   * denoted by uA_, lA_, uB_, and lB_, respectively
   */
   void Block::setupRadialLaplacian(Domain& domain)
   {
      // Preconditions
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(domainPtr_);
      UTIL_CHECK(domain.nr() > 0);
      UTIL_CHECK(domain.nz() > 0);

      int nr = domain.nr();
      double dr = domain.dr();
      double radius = domain.radius();

      // Check that work arrays are allocated, with correct sizes
      UTIL_CHECK(dA_.capacity() == nr);
      UTIL_CHECK(uA_.capacity() == nr - 1);
      UTIL_CHECK(lA_.capacity() == nr - 1);
      UTIL_CHECK(dB_.capacity() == nr);
      UTIL_CHECK(uB_.capacity() == nr - 1);
      UTIL_CHECK(lB_.capacity() == nr - 1);


      // Second derivative terms in matrix A
      double halfDs = 0.5*ds_;
      double db = kuhn()/dr;
      double c1 = halfDs*db*db/6.0;
      //double c2 = 2.0*c1;

      double halfDr = 0.5*dr;
      double x, rp, rm;

      // First row: x = xMin
      rp = 2.0*c1;
      dA_[0] += 2.0*rp;
      uA_[0] = -2.0*rp;

      // Interior rows
      for (int i = 1; i < nr - 1; ++i) {
         x = dr*i;
         rm = 1.0 - halfDr/x;
         rp = 1.0 + halfDr/x;
         rm *= c1;
         rp *= c1;
         dA_[i] += rm + rp;
         uA_[i] = -rp;
         lA_[i-1] = -rm;
      }

      // Last row: x = radius
      rm = 1.0 - halfDr/radius;
      rm *= c1;
      dA_[nr-1] += 2.0*rm;
      lA_[nr-2] = -2.0*rm;

      // Construct matrix B 
      for (int i = 0; i < nr; ++i) {
         dB_[i] = -dA_[i];
      }
      for (int i = 0; i < nr - 1; ++i) {
         uB_[i] = -uA_[i];
      }
      for (int i = 0; i < nr - 1; ++i) {
         lB_[i] = -lA_[i];
      }

      // Add diagonal identity terms to matrices A and B
      for (int i = 0; i < nr; ++i) {
         dA_[i] += 1.0;
         dB_[i] += 1.0;
      }

      // Compute the LU decomposition of matrix A 
      solver_.computeLU(dA_, uA_, lA_);
   }

   void Block::setupAxialLaplacian(Domain& domain) 
   {
      int nz = domain.nz();

   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Block::computeConcentration(double prefactor)
   {
      // Preconditions
      UTIL_CHECK(domain().nr() > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(dA_.capacity() == domain().nr());
      UTIL_CHECK(dB_.capacity() == domain().nr());
      UTIL_CHECK(uA_.capacity() == domain().nr() -1);
      UTIL_CHECK(uB_.capacity() == domain().nr() -1);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == domain().nr()) 

      // Initialize cField to zero at all points
      int i;
      int nr = domain().nr();
      for (i = 0; i < nr; ++i) {
         cField()[i] = 0.0;
      }

      Propagator const & p0 = propagator(0);
      Propagator const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (i = 0; i < nr; ++i) {
         cField()[i] += 0.5*p0.q(0)[i]*p1.q(ns_ - 1)[i];
      }
      for (int j = 1; j < ns_ - 1; ++j) {
         for (i = 0; i < nr; ++i) {
            cField()[i] += p0.q(j)[i]*p1.q(ns_ - 1 - j)[i];
         }
      }
      for (i = 0; i < nr; ++i) {
         cField()[i] += 0.5*p0.q(ns_ - 1)[i]*p1.q(0)[i];
      }

      // Normalize
      prefactor *= ds_;
      for (i = 0; i < nr; ++i) {
         cField()[i] *= prefactor;
      }

   }

   #if 0
   /*
   * Propagate solution by one step.
   *
   * This function implements one step of the Crank-Nicholson algorithm.
   * To do so, it solves A q(i+1) = B q(i), where A and B are constant 
   * matrices defined in the documentation of the setupStep() function.
   */
   void Block::step(const QField& q, QField& qNew)
   {
      int nr = domain().nr();
      v_[0] = dB_[0]*q[0] + uB_[0]*q[1];
      for (int i = 1; i < nr - 1; ++i) {
         v_[i] = dB_[i]*q[i] + lB_[i-1]*q[i-1] + uB_[i]*q[i+1];
      }
      v_[nr - 1] = dB_[nr-1]*q[nr-1] + lB_[nr-2]*q[nr-2];
      solver_.solve(v_, qNew);
   }
   #endif

}
}
