#ifndef PSCF_TRIDIAGONAL_SOLVER_H
#define PSCF_TRIDIAGONAL_SOLVER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

#include <string>
#include <iostream>


namespace Pscf 
{

   using namespace Util;

   /**
   * Solver for Ax=b with tridiagonal matrix A.
   *
   * \ingroup Pscf_Math_Module
   */  
   class TridiagonalSolver
   {
   public:

      /**
      * Constructor.
      */
      TridiagonalSolver();

      /**
      * Destructor.
      */
      ~TridiagonalSolver();

      /**
      * Allocate memory.
      *
      * \param n dimension of n x n square array.
      */
      void allocate(int n);

      /**
      * Compute LU decomposition of a symmetric tridiagonal matrix.
      *
      * \param d diagonal elements of n x n matrix matrix (0,..,n-1)
      * \param u upper off-diagonal elements (0,..,n-2)
      */
      void computeLU(const DArray<double>& d, const DArray<double>& u);

      /**
      * Compute LU decomposition of a general tridiagonal matrix.
      *
      * \param d diagonal elements of n x n matrix matrix (0,..,n-1)
      * \param u upper off-diagonal elements (0,..,n-2)
      * \param l lower off-diagonal elements (0,..,n-2)
      */
      void computeLU(const DArray<double>& d, 
                     const DArray<double>& u,
                     const DArray<double>& l);

      /**
      * Evaluate product Ab = x for known b to compute x.
      *
      * \param b known vector to be multiplied (input)
      * \param x result of multiplication Ab = x (output)
      */
      void multiply(const DArray<double>& b, DArray<double>& x);

      /**
      * Solve Ax = b for known b to compute x.
      *
      * \param b known vector on RHS (input)
      * \param x unknown solution vector of Ax = b (output)
      */
      void solve(const DArray<double>& b, DArray<double>& x);

   private:

      // Diagonal elements
      DArray<double> d_;

      // Upper off-diagonal elements (unmodified by computeLU)
      DArray<double> u_;

      // Lower off-diagonal elements (replaced by multipliers)
      DArray<double> l_;

      // Work space.
      DArray<double> y_;

      int n_;

      // Apply Gauss elimination to private arrays d_, u_, l_.
      void gaussElimination();

   };

}
#endif
