/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>
#include <utility>
#include <cmath>

#include "TridiagonalSolver.h"

namespace Pscf
{
  
   /*
   * Constructor.
   */
   TridiagonalSolver::TridiagonalSolver()
    : hasSwap_(false)
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
      u2_.allocate(n-1);
      l_.allocate(n-1);
      swap_.allocate(n-1);
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
      // Gauss elimination with threshold partial pivoting: rows i and
      // i+1 are swapped if the pivot d_[i] is small, |d_[i]| less than
      // 0.1*|l_[i]|, which adds an element u2_[i] in the second upper
      // off-diagonal of U. Without a swap, the multiplier is at most 10
      // in magnitude, and the result is that of elimination without
      // pivoting.
      double q, temp;
      hasSwap_ = false;
      for (int i = 0; i < n_ - 1; ++i) {
         if (i < n_ - 2) {
            u2_[i] = 0.0;
         }
         swap_[i] = (std::abs(d_[i]) < 0.1*std::abs(l_[i]));
         if (swap_[i]) {
            hasSwap_ = true;
            q = d_[i]/l_[i];
            d_[i] = l_[i];
            temp = u_[i];
            u_[i] = d_[i+1];
            d_[i+1] = temp - q*d_[i+1];
            if (i < n_ - 2) {
               u2_[i] = u_[i+1];
               u_[i+1] = -q*u_[i+1];
            }
         } else {
            q = l_[i]/d_[i];
            d_[i+1] -= q*u_[i];
         }
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
      // Without row swaps (the usual case)
      if (!hasSwap_) {
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
         return;
      }

      // Compute Ub = y
      for (int i = 0; i < n_ - 1; ++i) {
         y_[i] = d_[i]*b[i] + u_[i]*b[i+1];
         if (i < n_ - 2) {
            y_[i] += u2_[i]*b[i+2];
         }
      }
      y_[n_ - 1] = d_[n_ - 1]*b[n_ - 1];

      // Compute x = P_0 L_0 P_1 L_1 ... P_{n-2} L_{n-2} y, in which
      // L_i adds l_[i]*y_[i] to y_[i+1], and P_i swaps y_[i] and
      // y_[i+1] if rows i and i+1 were swapped
      for (int i = n_ - 2; i >= 0; --i) {
         y_[i+1] += l_[i]*y_[i];
         if (swap_[i]) {
            std::swap(y_[i], y_[i+1]);
         }
      }
      for (int i = 0; i < n_; ++i) {
         x[i] = y_[i];
      }
   }

   /*
   * Solve Ax = b for x, given known b.
   */
   void TridiagonalSolver::solve(const DArray<double>& b, DArray<double>& x)
   {
       // Without row swaps (the usual case)
       if (!hasSwap_) {
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
          return;
       }

       // Solve Ly = Pb by forward substitution, with the row swaps.
       for (int i = 0; i < n_; ++i) {
          y_[i] = b[i];
       }
       for (int i = 0; i < n_ - 1; ++i) {
          if (swap_[i]) {
             std::swap(y_[i], y_[i+1]);
          }
          y_[i+1] -= l_[i]*y_[i];
       } 

       // Solve Ux = y by back substitution.
       x[n_ - 1] = y_[n_ - 1]/d_[n_ - 1];
       for (int i = n_ - 2; i >= 0; --i) {
          x[i] = y_[i] - u_[i]*x[i+1];
          if (i < n_ - 2) {
             x[i] -= u2_[i]*x[i+2];
          }
          x[i] /= d_[i];
       }
   }

}
