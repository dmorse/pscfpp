/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

#include "TridiagonalSolver.h"

namespace Fd1d
{
  
   TridiagonalSolver::TridiagonalSolver()
   {}

   /*
   * Allocate memory.
   */
   void TridiagonalSolver::allocate(int n)
   {
      d_.allocate(n);
      u_.allocate(n-1);
      l_.allocate(n-1);
      y_.allocate(n);
      n_ = n;
   }

   /*
   * Compute the LU decomposition for later use.
   */
   void TridiagonalSolver::computeLU(const DArray<double>& d, const DArray<double>& u)
   {
      // Copy to local arrays
      for (int i = 0; i < n_ - 1; ++i) {
         d_[i] = d[i];
         u_[i] = u[i];
         l_[i] = u[i];
      }
      d_[n_ - 1] = d[n_ - 1];

      // Gauss elimination
      double q;
      for (int i = 0; i < n_ - 1; ++i) {
         q = l_[i]/d_[i];
         d_[i+1] -= q*u_[i];
         l_[i] = q;
      }
   }

   /*
   * Solve Ax = b.
   */
   void TridiagonalSolver::solve(const DArray<double>& b, DArray<double>& x)
   {
       // Solve Ly = b by forward substitution.
       y_[0] = b[0];
       for (int i = 1; i < n_; ++i) {
          y_[i] = b[i] - l_[i-1]*y_[i-1]; 
       } 

       // Solve Ux = y by back substitution.
       x[n_ - 1] = y_[n_ - 1];
       for (int i = n_ - 2; i >= 0; --i) {
          x[i] = (y_[i] - u_[i]*x[i+1])/d_[i]; 
       }
   }

}
