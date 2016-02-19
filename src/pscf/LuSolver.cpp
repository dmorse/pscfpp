/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

// #include <string>
// #include <iostream>

#include "LuSolver.h"
#include <gsl/gsl_linalg.h>

namespace Pscf
{
  
   LuSolver::LuSolver()
    : luPtr_(0),
      permPtr_(0),
      n_(0)
   {}

   /*
   * Allocate memory.
   */
   void LuSolver::allocate(int n)
   {
      luPtr_ = gsl_matrix_alloc(n, n);
      permPtr_ = gsl_permutation_alloc(n);
      n_ = n;
   }

   /*
   * Compute the LU decomposition for later use.
   */
   void LuSolver::computeLU(const Matrix<double>& A)
   {
      int i, j, k;
      k = 0;
      for (i = 0; i < n_; ++i) {
         for (j = 0; j < n_; ++i) {
            luPtr_->data[k] = A(i, j); 
         }
      }
      gsl_linalg_LU_decomp(luPtr_, permPtr_, &signum_);
   }

   /*
   * Solve Ax = b.
   */
   void LuSolver::solve(Array<double>& b, Array<double>& x)
   {
      // Associate gsl_vector b_ with Array b
      b_.size = n_;
      b_.stride = 1;
      b_.data = b.cArray();
      b_.block= 0;
      b_.owner= 0;

      // Associate gsl_vector x_ with Array x
      x_.size = n_;
      x_.stride = 1;
      x_.data = x.cArray();
      x_.block= 0;
      x_.owner= 0;

      gsl_linalg_LU_solve(luPtr_, permPtr_, &b_, &x_);
   }

}
