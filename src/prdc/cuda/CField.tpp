#ifndef PRDC_CUDA_C_FIELD_TPP
#define PRDC_CUDA_C_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/CField.h>

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;
   using namespace Pscf;

   /**
   * Default constructor.
   */
   template <int D>
   CField<D>::CField()
    : Field<cudaComplex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   CField<D>::~CField()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <int D>
   CField<D>::CField(const CField<D>& other)
    : Field<cudaComplex>(other),
      meshDimensions_(0)
   {
      meshDimensions_ = other.meshDimensions_;
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other Field is not allocated.
   * \throw Exception if both Fields are allocated with unequal capacities.
   *
   * \param other the rhs Field
   */
   template <int D>
   CField<D>& CField<D>::operator = (const CField<D>& other)
   {
      Field<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
}
}
#endif
