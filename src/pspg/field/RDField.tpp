#ifndef PSPG_R_DFIELD_TPP
#define PSPG_R_DFIELD_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RDField.h"

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RDField<D>::RDField()
    : DField<cudaReal>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RDField<D>::~RDField()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <int D>
   RDField<D>::RDField(const RDField<D>& other)
    : DField<cudaReal>(other),
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
   RDField<D>& RDField<D>::operator = (const RDField<D>& other)
   {
      DField<cudaReal>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
}
#endif
