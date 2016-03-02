/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LuSolver.h"
#include <gsl/gsl_linalg.h>

namespace Pscf
{
  
   LuSolver::LuSolver()
    : luPtr_(0),
      permPtr_(0),
      n_(0)
   {
      // Initialize gs_vector b_ 
      b_.size = 0;
      b_.stride = 1;
      b_.data = 0;
      b_.block = 0;
      b_.owner = 0;

      // Initialize gsl_vector x_ 
      x_.size = 0;
      x_.stride = 1;
      x_.data = 0;
      x_.block = 0;
      x_.owner = 0;
   }

   LuSolver::~LuSolver()
   {
      if (n_ > 0)  {
         if (luPtr_) gsl_matrix_free(luPtr_);
         if (permPtr_) gsl_permutation_free(permPtr_);
      }
   }

   /*
   * Allocate memory and set n.
   */
   void LuSolver::allocate(int n)
   {
      UTIL_CHECK(n > 0);
      UTIL_CHECK(n_ == 0);
      UTIL_CHECK(luPtr_ == 0);
      UTIL_CHECK(permPtr_ == 0);
      luPtr_ = gsl_matrix_alloc(n, n);
      permPtr_ = gsl_permutation_alloc(n);
      b_.size = n;
      x_.size = n;
      n_ = n;
   }

   /*
   * Compute the LU decomposition for later use.
   */
   void LuSolver::computeLU(const Matrix<double>& A)
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(A.capacity1() == n_);
      UTIL_CHECK(A.capacity2() == n_);

      int i, j;
      int k = 0;
      for (i = 0; i < n_;  ++i) {
         for (j = 0; j < n_; ++j) {
            luPtr_->data[k] = A(i,j);
            ++k;
         }
      }
      gsl_linalg_LU_decomp(luPtr_, permPtr_, &signum_);
   }

   /*
   * Solve Ax = b.
   */
   void LuSolver::solve(Array<double>& b, Array<double>& x)
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(b.capacity() == n_);
      UTIL_CHECK(x.capacity() == n_);

      // Associate gsl_vectors b_ and x_ with Arrays b and x
      b_.data = b.cArray();
      x_.data = x.cArray();

      // Solve system of equations
      gsl_linalg_LU_solve(luPtr_, permPtr_, &b_, &x_);

      // Destroy temporary associations
      b_.data = 0;
      x_.data = 0;
   }

}
