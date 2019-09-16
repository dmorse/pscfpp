#ifndef PSCF_SPACE_SYMMETRY_TPP
#define PSCF_SPACE_SYMMETRY_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpaceSymmetry.h"

namespace Pscf
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   SpaceSymmetry<D>::SpaceSymmetry()
    : R_(),
      t_()
   {
      int i, j;
      for (i = 0; i < D; ++i) {
         t_[i] = 0;
         for (j = 0; j < D; ++j) {
            if (i == j) {
               R_(i, j) = 1;
            } else {
               R_(i, j) = 0;
            }
         }
      }
   }

   /*
   * Copy constructor.
   */
   template <int D>
   SpaceSymmetry<D>::SpaceSymmetry(const SpaceSymmetry<D>& other)
   {
      int i, j;
      for (i = 0; i < D; ++i) {
         t_[i] = other.t_[i];
         for (j = 0; j < D; ++j) {
            R_(i, j) = other.R_(i,j);
         }
      }
      normalize();
   }

   /*
   * Assignment operator.
   */
   template <int D>
   SpaceSymmetry<D>& 
   SpaceSymmetry<D>::operator = (const SpaceSymmetry<D>& other)
   {
      if (this != &other) {
         int i, j;
         for (i = 0; i < D; ++i) {
            t_[i] = other.t_[i];
            for (j = 0; j < D; ++j) {
               R_(i, j) = other.R_(i,j);
            }
         }
      }
      normalize();
      return *this;
   }

   /*
   * Shift translation to lie in range [0,1). 
   */
   template <int D>
   void SpaceSymmetry<D>::normalize()
   {
      for (int i = 0; i < D; ++i) {
         int num = t_[i].num();
         int den = t_[i].den();
         UTIL_ASSERT(den > 0);
         if (den == 1) {
            num = 0;
         } else {
            while (num < 0) {
               num += den;
            }
            while (num >= den) {
              num -= den;
            }
         }
         if (num != t_[i].num()) {
            t_[i] = Rational(num, den);
         }
      }
   }

   /*
   * Make identity (private static method).
   */
   template <int D>
   void SpaceSymmetry<D>::makeIdentity()
   {
      int i, j;
      for (i = 0; i < D; ++i) {
         identity_.t_[i] = 0;
         for (j = 0; j < D; ++j) {
            if (i == j) {
               identity_.R_(i, j) = 1;
            } else {
               identity_.R_(i, j) = 0;
            }
         }
      }
   }

   /*
   * Return inverse of this SpaceSymmetry<D>.
   */
   template <int D>
   SpaceSymmetry<D> SpaceSymmetry<D>::inverse() const
   {
      SpaceSymmetry<D> C;

      // Compute inverse of rotation matrix
      C.R_ = inverseRotation();

      // Compute translation -R^{-1}t
      int i, j;
      for (i = 0; i < D; ++i) {
         C.t_[i] = 0;
         for (j = 0; j < D; ++j) {
            C.t_[i] -= C.R_(i, j)*t_[j];
         }
      }
      C.normalize();

      //UTIL_CHECK(C.determinant()*determinant() == 1);

      return C;
   }

}
#endif
