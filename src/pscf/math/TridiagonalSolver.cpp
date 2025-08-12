/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

#include "TridiagonalSolver.h"

namespace Pscf
{
  
   /*
   * Constructor.
   */
   TridiagonalSolver::TridiagonalSolver()
   {}

   /*
   * Destructor.
   */
   TridiagonalSolver::~TridiagonalSolver()
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
   * Compute the LU decomposition of a symmetric matrix for later use.
   */
   void TridiagonalSolver::computeLU(const DArray<double>& d, 
                                     const DArray<double>& u)
   {
      // Copy to local arrays
      for (int i = 0; i < n_ - 1; ++i) {
         d_[i] = d[i];
         u_[i] = u[i];
         l_[i] = u[i];
      }
      d_[n_ - 1] = d[n_ - 1];
      gaussElimination();
   }

   /*
   * Compute the LU decomposition of a symmetric matrix for later use.
   */
   void TridiagonalSolver::computeLU(const DArray<double>& d, 
                                     const DArray<double>& u,
                                     const DArray<double>& l)
   {
      // Copy to local arrays
      for (int i = 0; i < n_ - 1; ++i) {
         d_[i] = d[i];
         u_[i] = u[i];
         l_[i] = l[i];
      }
      d_[n_ - 1] = d[n_ - 1];
      gaussElimination();
   }

   void TridiagonalSolver::gaussElimination()
   {
      // Gauss elimination
      double q;
      for (int i = 0; i < n_ - 1; ++i) {
         q = l_[i]/d_[i];
         d_[i+1] -= q*u_[i];
         l_[i] = q;
      }

      #if 0
      std::cout << "\n";
      for (int i=0; i < n_; ++i) {
         std::cout << "  " << d_[i];
      }
      std::cout << "\n";
      for (int i=0; i < n_ - 1; ++i) {
         std::cout << "  " << u_[i];
      }
      std::cout << "\n";
      for (int i=0; i < n_ - 1; ++i) {
         std::cout << "  " << l_[i];
      }
      std::cout << "\n";
      #endif
   }

   /*
   * Compute matrix product Ab = x for x, given known vector b.
   */
   void TridiagonalSolver::multiply(const DArray<double>& b, DArray<double>& x)
   {
      // Compute Ub = y
      for (int i = 0; i < n_ - 1; ++i) {
         y_[i] = d_[i]*b[i] + u_[i]*b[i+1];
      }
      y_[n_ - 1] = d_[n_ - 1]*b[n_ - 1];

      // Compute Ly = x
      x[0] = y_[0];
      for (int i = 1; i < n_; ++i) {
         x[i] = y_[i] + l_[i-1]*y_[i-1];
      }
   }

   /*
   * Solve Ax = b for x, given known b.
   */
   void TridiagonalSolver::solve(const DArray<double>& b, DArray<double>& x)
   {
       // Solve Ly = b by forward substitution.
       y_[0] = b[0];
       for (int i = 1; i < n_; ++i) {
          y_[i] = b[i] - l_[i-1]*y_[i-1]; 
       } 

       // Solve Ux = y by back substitution.
       x[n_ - 1] = y_[n_ - 1]/d_[n_ - 1];
       for (int i = n_ - 2; i >= 0; --i) {
          x[i] = (y_[i] - u_[i]*x[i+1])/d_[i]; 
       }
   }

}
