#ifndef PSCF_LU_SOLVER_H
#define PSCF_LU_SOLVER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
      * Compute inverse of matrix A.
      *
      * \param inv inverse of matrix A (output)
      */  
      void inverse (Matrix<double>& inv); 
 
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
