#ifndef PSCF_LU_SOLVER_H
#define PSCF_LU_SOLVER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/containers/Matrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>

namespace Pscf 
{

   using namespace Util;

   /**
   * Solve Ax=b by LU decomposition of A.
   *
   * This class is a simple wrapper for the functions provided by
   * the Gnu Scientific Library (GSL).
   *
   * \ingroup Pscf_Math_Module
   */  
   class LuSolver
   {
   public:

      /**
      * Constructor.
      */
      LuSolver();

      /**
      * Destructor.
      */
      ~LuSolver();

      /**
      * Allocate memory.
      *
      * \param n dimension of n x n square array.
      */
      void allocate(int n);

      /**
      * Compute the LU decomposition for later use.
      *
      * The decomposition is completed also for a singular matrix A,
      * with a zero pivot, which isSingular() then reports.
      *
      * \param A the square matrix A in problem Ax=b.
      */
      void computeLU(const Matrix<double>& A);

      /**
      * Solve Ax = b for known b to compute x.
      *
      * \param b the RHS vector
      * \param x the solution vector
      */
      void solve(Array<double>& b, Array<double>& x);

      /**
      * Solve Ax = b by a basic solution, also for a singular matrix A.
      *
      * The component of x of each zero pivot of the LU decomposition
      * is set to zero, and the equation of that pivot is not used. The
      * result solves Ax = b if b lies in the range of A and the number
      * of zero pivots equals the rank deficiency of A, i.e., if each
      * zero pivot belongs to a column that depends on the preceding
      * columns (as for a matrix with a single null direction).
      * Otherwise it need not solve the system, so callers must check
      * the residual. Use solve() for a matrix that is not singular.
      *
      * \param b the RHS vector
      * \param x the solution vector (must not be b)
      */
      void solveBasic(Array<double> const & b, Array<double>& x) const;

      /** 
      * Compute inverse of matrix A.
      *
      * \param inv inverse of matrix A (output)
      */  
      void inverse (Matrix<double>& inv); 

      /**
      * Does the last LU decomposition have a zero pivot?
      *
      * If so, solve() and inverse() abort through the GSL error handler.
      */
      bool isSingular() const;
 
   private:

      /// RHS vector of Ax=b.
      gsl_vector b_;

      /// Solution vector of Ax=b.
      gsl_vector x_;

      /// Pointer to LU decomposition matrix.
      gsl_matrix* luPtr_;

      /// Pointer to inverse matrix.
      gsl_matrix* gMatInverse_;

      /// Pointer to permutation in LU decomposition.
      gsl_permutation* permPtr_;

      /// Sign of permutation in LU decomposition.
      int signum_;

      /// Number of rows and columns in matrix.
      int n_;

   };

}
#endif
