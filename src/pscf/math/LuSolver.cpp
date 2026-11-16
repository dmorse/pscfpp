/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "LuSolver.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_version.h>
#include <cmath>
#include <vector>

// Since GSL 2.6, gsl_linalg_LU_decomp is recursive (see computeLU)
#if defined(GSL_MAJOR_VERSION) && defined(GSL_MINOR_VERSION)
#if GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 6)
#define PSCF_LU_RECURSIVE
#endif
#endif

namespace Pscf
{

#ifdef PSCF_LU_RECURSIVE
   namespace {

      /*
      * LU decomposition with partial pivoting of an M x N matrix A,
      * M >= N, P A = L U, where rows j and ipiv[j] were swapped in
      * step j (as in LAPACK getrf).
      *
      * This is the recursive algorithm of gsl_linalg_LU_decomp of
      * GSL 2.8, from the functions LU_decomp_L3 and LU_decomp_L2 of
      * the GSL file linalg/lu.c,
      *
      *   Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009
      *   Gerard Jungman, Brian Gough
      *   Copyright (C) 2019, 2021 Patrick Alken
      *
      * distributed under the GNU General Public License, version 3 or
      * later, like PSCF. It does the same floating point operations,
      * except that the decomposition is completed past a zero pivot:
      * GSL 2.8 returns from each level of the recursion as soon as a
      * block reports a zero pivot, leaving the later columns
      * unfactored and their row swaps unset or unshifted (GSL 2.6 and
      * 2.7 divide by the zero pivot, which fills the remaining
      * submatrix with NaN).
      */
      void decomposeLU(gsl_matrix* A, size_t* ipiv)
      {
         const size_t M = A->size1;
         const size_t N = A->size2;
         size_t i, j;
         if (N <= 24) {
            // Unblocked algorithm (LAPACK getf2)
            for (j = 0; j < N; ++j) {
               gsl_vector_view c = gsl_matrix_subcolumn(A, j, j, M - j);
               ipiv[j] = j + gsl_blas_idamax(&c.vector);
               if (gsl_matrix_get(A, ipiv[j], j) != 0.0) {
                  if (ipiv[j] != j) gsl_matrix_swap_rows(A, j, ipiv[j]);
                  double ajj = gsl_matrix_get(A, j, j);
                  if (j < M - 1 && std::fabs(ajj) >= GSL_DBL_MIN) {
                     c = gsl_matrix_subcolumn(A, j, j + 1, M - j - 1);
                     gsl_blas_dscal(1.0 / ajj, &c.vector);
                  } else {
                     for (i = j + 1; i < M; ++i) {
                        *gsl_matrix_ptr(A, i, j) /= ajj;
                     }
                  }
               }
               if (j < N - 1) {
                  gsl_vector_view l, u;
                  l = gsl_matrix_subcolumn(A, j, j + 1, M - j - 1);
                  u = gsl_matrix_subrow(A, j, j + 1, N - j - 1);
                  gsl_matrix_view A22;
                  A22 = gsl_matrix_submatrix(A, j+1, j+1, M-j-1, N-j-1);
                  gsl_blas_dger(-1.0, &l.vector, &u.vector, &A22.matrix);
               }
            }
            return;
         }

         // Recursion: columns 0 to N1 - 1, then the submatrix A22
         const size_t N1 = ((N + 8) / 16) * 8;
         const size_t N2 = N - N1;
         gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
         gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, N1, N2);
         gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, M - N1, N1);
         gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, M - N1, N2);
         gsl_matrix_view AL = gsl_matrix_submatrix(A, 0, 0, M, N1);
         gsl_matrix_view AR = gsl_matrix_submatrix(A, 0, N1, M, N2);
         decomposeLU(&AL.matrix, ipiv);
         for (i = 0; i < N1; ++i) {
            if (ipiv[i] != i) gsl_matrix_swap_rows(&AR.matrix, i, ipiv[i]);
         }
         gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                        1.0, &A11.matrix, &A12.matrix);
         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A21.matrix,
                        &A12.matrix, 1.0, &A22.matrix);
         decomposeLU(&A22.matrix, ipiv + N1);
         for (i = N1; i < N; ++i) {
            if (ipiv[i] != i - N1) {
               gsl_matrix_swap_rows(&A21.matrix, i - N1, ipiv[i]);
            }
            ipiv[i] += N1;
         }
      }

   }
#endif
  
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
         if (gMatInverse_) gsl_matrix_free(gMatInverse_);
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
      gMatInverse_ = gsl_matrix_alloc(n, n);
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

      // LU decomposition with partial pivoting, P A = L U, in the
      // format of gsl_linalg_LU_decomp, completed also for a singular
      // matrix, with a zero pivot that isSingular() reports. Since GSL
      // 2.6, gsl_linalg_LU_decomp decomposes matrices with more than
      // 24 columns recursively and does not complete the decomposition
      // of a singular matrix, so its algorithm is used here (function
      // decomposeLU above); for a nonsingular matrix, the result is that
      // of gsl_linalg_LU_decomp, bit for bit. GSL 2.5 and earlier use
      // an unblocked algorithm, which continues past a zero pivot.
#ifdef PSCF_LU_RECURSIVE
      std::vector<size_t> ipiv(n_);
      decomposeLU(luPtr_, ipiv.data());
      gsl_permutation_init(permPtr_);
      signum_ = 1;
      for (i = 0; i < n_; ++i) {
         if (ipiv[i] != (size_t) i) {
            gsl_permutation_swap(permPtr_, i, ipiv[i]);
            signum_ = -signum_;
         }
      }
#else
      gsl_linalg_LU_decomp(luPtr_, permPtr_, &signum_);
#endif
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

   /*
   * Does the LU decomposition have a zero pivot (GSL's singular test)?
   */
   bool LuSolver::isSingular() const
   {
      UTIL_CHECK(n_ > 0);
      for (int i = 0; i < n_; ++i) {
         if (gsl_matrix_get(luPtr_, i, i) == 0.0) return true;
      }
      return false;
   }

   /*
   * Solve Ax = b by a basic solution, with zero pivots allowed.
   */
   void LuSolver::solveBasic(Array<double> const & b, Array<double>& x)
   const
   {
      UTIL_CHECK(n_ > 0);
      UTIL_CHECK(b.capacity() == n_);
      UTIL_CHECK(x.capacity() == n_);
      UTIL_CHECK(b.cArray() != x.cArray());

      int i, j;
      double sum, pivot;

      // Row permutation, x = P b
      for (i = 0; i < n_; ++i) {
         x[i] = b[gsl_permutation_get(permPtr_, i)];
      }

      // Forward substitution, L y = P b (L has a unit diagonal)
      for (i = 1; i < n_; ++i) {
         sum = x[i];
         for (j = 0; j < i; ++j) {
            sum -= gsl_matrix_get(luPtr_, i, j)*x[j];
         }
         x[i] = sum;
      }

      // Back substitution, U x = y, with x[i] = 0 for a zero pivot
      for (i = n_ - 1; i >= 0; --i) {
         pivot = gsl_matrix_get(luPtr_, i, i);
         if (pivot == 0.0) {
            x[i] = 0.0;
         } else {
            sum = x[i];
            for (j = i + 1; j < n_; ++j) {
               sum -= gsl_matrix_get(luPtr_, i, j)*x[j];
            }
            x[i] = sum/pivot;
         }
      }
   }

   /*  
   * Find inverse 
   */  
   void LuSolver::inverse(Matrix<double>& inv)
   {   
      UTIL_CHECK(n_ > 0); 

      gMatInverse_->data = inv.cArray();
      gsl_linalg_LU_invert(luPtr_, permPtr_, gMatInverse_);

   } 

}
