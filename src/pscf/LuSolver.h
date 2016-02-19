#ifndef PSCF_LU_SOLVER_H
#define PSCF_LU_SOLVER_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
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
   * \ingroup Pscf_Base_Module
   */  
   class LuSolver
   {
   public:

      /**
      * Constructor.
      */
      LuSolver();

      /**
      * Allocate memory.
      *
      * \param n dimension of n x n square array.
      */
      void allocate(int n);

      /**
      * Compute the LU decomposition for later use.
      */
      void computeLU(const Matrix<double>& A);

      /**
      * Evaluate product Ab = x for known b to compute x.
      */
      void multiply(const Array<double>& b, Array<double>& x);

      /**
      * Solve Ax = b for known b to compute x.
      */
      void solve(Array<double>& b, Array<double>& x);

   private:

      // RHS b vector
      gsl_vector b_;

      // Solution vector
      gsl_vector x_;

      // Pointer to LU decomposition
      gsl_matrix* luPtr_;

      // Pointer to permutation in LU decomposition
      gsl_permutation* permPtr_;

      /// Pointer of permutation in LU decomposition
      int signum_;

      // Number of rows and columns in matrix
      int n_;

   };

}
#endif
