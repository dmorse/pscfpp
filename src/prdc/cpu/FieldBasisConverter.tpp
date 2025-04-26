#ifndef PRDC_CPU_FIELD_BASIS_CONVERTER_TPP
#define PRDC_CPU_FIELD_BASIS_CONVERTER_TPP

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldBasisConverter.h"
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FieldBasisConverter<D>::FieldBasisConverter()
    : basis_(),
      nMonomer_(0)
   {}
      
   /*
   * Constructor, creates basis.
   */
   template <int D>
   FieldBasisConverter<D>::FieldBasisConverter(DMatrix<double> basis)
    : basis_(),
      nMonomer_(0)
   {  setBasis(basis); }

   /*
   * Destructor.
   */
   template <int D>
   FieldBasisConverter<D>::~FieldBasisConverter()
   {}

   /*
   * Set the basis.
   */
   template <int D>
   void FieldBasisConverter<D>::setBasis(DMatrix<double> basis)
   {
      UTIL_CHECK(basis.isAllocated());
      UTIL_CHECK(basis.capacity1() > 1);
      UTIL_CHECK(basis.capacity1() == basis.capacity2());

      basis_ = basis;
      nMonomer_ = basis.capacity1();

      UTIL_CHECK(basis_.capacity1() == nMonomer_);
      UTIL_CHECK(basis_.capacity2() == nMonomer_);
   }

   template <int D> 
   double FieldBasisConverter<D>::maxBasisError(double normSq) const
   {
      UTIL_CHECK(nMonomer_ > 1);
      UTIL_CHECK(normSq > 0.0);

      double error = 0.0;
      double maxError = 0.0;
      int i, j, k;
      for (i = 0; i < nMonomer_; ++i) {
         for (j = 0; j <= i; ++j) {
            error = 0.0;
            for (k = 0; k < nMonomer_; ++k) {
               error +=  basis_(i,k)*basis_(j,k);
            }
            if (i == j) {
               error = error - normSq;
            }
            error = std::abs(error);
            if (error > maxError) {
               maxError = error;
            }
         }
      }
      return maxError;
   }
   
   /*
   * Compute pointwise components of r-grid fields in basis.
   */
   template <int D> 
   void 
   FieldBasisConverter<D>::convertToBasis(DArray< RField<D> > const & in,
                                          DArray< RField<D> > & out,
                                          double prefactor) const
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 1);
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() == nMonomer_);
      UTIL_CHECK(out.capacity() == nMonomer_);
      int meshSize = in[0].capacity();
      for (int i=0; i < nMonomer_; ++i) {
         UTIL_CHECK(in[i].capacity() == meshSize);
         UTIL_CHECK(out[i].capacity() == meshSize);
      }

      double vec;
      int i, j, k;

      // Loop over components in basis
      for (i = 0; i < nMonomer_; ++i) {
         RField<D>& wo = out[i];

         // Initialize wo = out[i] to zero
         for (k = 0; k < meshSize; ++k) {
            wo[k] = 0.0; 
         }

         // Loop over monomer types and mesh points
         for (j = 0; j < nMonomer_; ++j) {
            vec = basis_(i, j)*prefactor;
            RField<D> const & wi = in[j];
            for (k = 0; k < meshSize; ++k) {
               wo[k] += vec*wi[k];
            }
         } // for j < nMonomer_

      } // for i < nMonomer_

   }

   /*
   * Compute r-grid fields from components in basis.
   */
   template <int D> 
   void 
   FieldBasisConverter<D>::convertFromBasis(DArray<RField<D>> const & in,
                                            DArray<RField<D>> & out,
                                            double prefactor) const
   {
      // Preconditions
      UTIL_CHECK(nMonomer_ > 1);
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() == nMonomer_);
      UTIL_CHECK(out.capacity() == nMonomer_);
      int meshSize = in[0].capacity();
      for (int i=0; i < nMonomer_; ++i) {
         UTIL_CHECK(in[i].capacity() == meshSize);
         UTIL_CHECK(out[i].capacity() == meshSize);
      }

      double vec;
      int i, j, k;

      // Loop over monomer types
      for (i = 0; i < nMonomer_; ++i) {
         RField<D>& wo = out[i];

         // Initialize wo = out[i] to zero
         for (k = 0; k < meshSize; ++k) {
            wo[k] = 0.0; 
         }

         // Loop over components in basis 
         for (j = 0; j < nMonomer_; ++j) {
            vec = basis_(j, i)*prefactor;
            RField<D> const & wi = in[j];
            // Loop over mesh points
            for (k = 0; k < meshSize; ++k) {
               wo[k] += vec*wi[k];
            }
         } // for j < nMonomer_

      } // for i < nMonomer_

   }

}
}
}
#endif
