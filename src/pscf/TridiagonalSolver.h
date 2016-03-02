#ifndef PSCF_TRIDIAGONAL_SOLVER_H
#define PSCF_TRIDIAGONAL_SOLVER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

#include <string>
#include <iostream>


namespace Pscf 
{

   using namespace Util;

   /**
   * Solver for Ax=b with tridiagonal, symmetric matrix A.
   *
   * \ingroup Pscf_Base_Module
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
      * Compute the LU decomposition for later use.
      */
      void computeLU(const DArray<double>& d, const DArray<double>& u);

      /**
      * Evaluate product Ab = x for known b to compute x.
      */
      void multiply(const DArray<double>& b, DArray<double>& x);

      /**
      * Solve Ax = b for known b to compute x.
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
   };

}
#endif
