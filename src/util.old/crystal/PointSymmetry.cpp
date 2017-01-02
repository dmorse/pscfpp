

#include "PointSymmetry.h"

#include <util/format/Int.h>
#include <util/space/IntVector.h>

namespace Util 
{

   // Define static members

   PointSymmetry PointSymmetry::identity_;
   bool          PointSymmetry::hasIdentity_ = false;

   // Define methods

   /*
   * Default constructor.
   */
   PointSymmetry::PointSymmetry()
    : R_()
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            R_(i, j) = 0;
         }
      }
   }

   /*
   * Copy constructor.
   */
   PointSymmetry::PointSymmetry(const PointSymmetry& other)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            R_(i, j) = other.R_(i,j);
         }
      }
   }

   /*
   * Make identity (private static method).
   */
   void PointSymmetry::makeIdentity()
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            if (i == j) {
               identity_.R_(i, j) = 1;
            } else {
               identity_.R_(i, j) = 0;
            }
         }
      }
   }

   /*
   * Return inverse of this PointSymmetry.
   */
   PointSymmetry PointSymmetry::inverse() const
   {
      PointSymmetry C;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            C.R_(i, j) = R_(j, i);
         }
      }
      return C;
   }

   /*
   * Assignment operator.
   */
   PointSymmetry& PointSymmetry::operator = (const PointSymmetry& other)
   {
      if (this != &other) {
         int i, j;
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < Dimension; ++j) {
               R_(i, j) = other.R_(i,j);
            }
         }
      }
      return *this;
   }

   /*
   * Equality operator.
   */
   bool operator == (const PointSymmetry& A, const PointSymmetry& B)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            if (A.R_(i, j) != B.R_(i,j)) {
               return false;
            }
         }
      }
      return true;
   }

   /*
   * Group multipication operator for PointSymmetry objects.
   */
   PointSymmetry operator * (const PointSymmetry& A, const PointSymmetry& B)
   {
      PointSymmetry C;
      int i, j, k;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            C.R_(i, j) = 0;
            for (k = 0; k < Dimension; ++k) {
               C.R_(i, j) += A.R_(i, k)*B.R_(k, j);
            }
         }
      }
      return C;
   }

   /*
   * Matrix / IntVector multiplication.
   */
   IntVector operator * (const PointSymmetry& S, const IntVector& V)
   {
      IntVector U;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         U[i] = 0;
         for (j = 0; j < Dimension; ++j) {
            U[i] += S.R_(i,j)*V[j];
         }
      }
      return U;
   }

   /*
   * IntVector / Matrix multiplication.
   */
   IntVector operator * (const IntVector& V, const PointSymmetry& S)
   {
      IntVector U;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         U[i] = 0;
         for (j = 0; j < Dimension; ++j) {
            U[i] += V[j]*S.R_(j,i);
         }
      }
      return U;
   }

   /*
   * Output stream inserter for a PointSymmetry.
   */ 
   std::ostream& operator << (std::ostream& out, const PointSymmetry& A)
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            out << Int(A.R_(i,j));
         }
         out << std::endl;
      }
      return out;
   }

}
