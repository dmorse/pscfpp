/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
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
   {}

   LuSolver::~LuSolver()
   {
      if (n_ > 0)  {
         if (luPtr_) gsl_matrix_free(luPtr_);
         if (permPtr_) gsl_permutation_free(permPtr_);
      }
   }

   /*
   * Allocate memory.
   */
   void LuSolver::allocate(int n)
   {
      UTIL_CHECK(n > 0);
      luPtr_ = gsl_matrix_alloc(n, n);
      permPtr_ = gsl_permutation_alloc(n);
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

      // Associate gsl_vector b_ with Array b
      b_.size = n_;
      b_.stride = 1;
      b_.data = b.cArray();
      b_.block = 0;
      b_.owner = 0;

      // Associate gsl_vector x_ with Array x
      x_.size = n_;
      x_.stride = 1;
      x_.data = x.cArray();
      x_.block = 0;
      x_.owner = 0;

      gsl_linalg_LU_solve(luPtr_, permPtr_, &b_, &x_);
   }

}
