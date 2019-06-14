#include "SpaceSymmetry.h"

#include <util/format/Int.h>

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
      return *this;
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

   #if 0
   /*
   * Return inverse of rotation matrix.
   *
   * The general template is unimplemented, and thus throws and error.
   */
   template <int D>
   SpaceSymmetry<D> SpaceSymmetry<D>::inverseRotation() const
   {
      UTIL_THROW("Unimplemented SpaceSymmetry<D>::inverseRotation()");
      SpaceSymmetry<D> C;
      return C;
   }
   #endif

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

      return C;
   }

   /*
   * Equality operator.
   */
   template <int D>
   bool operator == (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B)
   {
      int i, j;
      for (i = 0; i < D; ++i) {
         if (A.t_[i] != B.t_[i]) {
            return false;
         }
         for (j = 0; j < D; ++j) {
            if (A.R_(i, j) != B.R_(i,j)) {
               return false;
            }
         }
      }
      return true;
   }

   /*
   * Group multipication operator for SpaceSymmetry<D> objects.
   */
   template <int D>
   SpaceSymmetry<D> 
   operator * (const SpaceSymmetry<D>& A, const SpaceSymmetry<D>& B)
   {
      SpaceSymmetry<D> C;
      int i, j, k;

      // Compute rotation matrix (matrix product)
      for (i = 0; i < D; ++i) {
         for (j = 0; j < D; ++j) {
            C.R_(i, j) = 0;
            for (k = 0; k < D; ++k) {
               C.R_(i, j) += A.R_(i, k)*B.R_(k, j);
            }
         }
      }

      // Compute translation vector
      for (i = 0; i < D; ++i) {
         C.t_[i] = A.t_[i];
      }
      for (i = 0; i < D; ++i) {
         for (j = 0; j < D; ++j) {
            C.t_[i] += A.R_(i, j)*B.t_[j];
         }
      }

      return C;
   }

   /*
   * Matrix / IntVec<D> multiplication.
   */
   template <int D>
   IntVec<D> operator * (const SpaceSymmetry<D>& S, const IntVec<D>& V)
   {
      IntVec<D> U;
      int i, j;
      for (i = 0; i < D; ++i) {
         U[i] = 0;
         for (j = 0; j < D; ++j) {
            U[i] += S.R_(i,j)*V[j];
         }
      }
      return U;
   }

   /*
   * IntVec<D> / Matrix multiplication.
   */
   template <int D>
   IntVec<D> operator * (const IntVec<D>& V, const SpaceSymmetry<D>& S)
   {
      IntVec<D> U;
      int i, j;
      for (i = 0; i < D; ++i) {
         U[i] = 0;
         for (j = 0; j < D; ++j) {
            U[i] += V[j]*S.R_(j,i);
         }
      }
      return U;
   }

   /*
   * Output stream inserter for a SpaceSymmetry<D>.
   */ 
   template <int D>
   std::ostream& operator << (std::ostream& out, const SpaceSymmetry<D>& A)
   {
      int i, j;
      for (i = 0; i < D; ++i) {
         for (j = 0; j < D; ++j) {
            out << " " << Int(A.R_(i,j),2);
         }
         out << std::endl;
      }
      for (i = 0; i < D; ++i) {
         out << "  " << A.t_[i];
      }
      out << std::endl;
      return out;
   }

}
