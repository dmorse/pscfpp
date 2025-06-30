#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"              // parent class template parameter
#include <prdc/cuda/RField.h>     // parent class template parameter
#include <prdc/field/MaskReal.h>  // parent class

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * A field to which the total monomer concentration is constrained.
   *
   * Please refer to the documentation of the Prdc::MaskReal base class
   * template for more complete API documentation for this class template.
   * The public interface of Rpg::Mask is identical to that of the base
   * class template Prdc::MaskReal. 
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Mask 
     : public Prdc::MaskReal< D, FieldIo<D>, Prdc::Cuda::RField<D> >
   {

   public:

      /// Base class typedef
      typedef Prdc::MaskReal< D, FieldIo<D>, Prdc::Cuda::RField<D> > Base;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocateBasis;
      using Base::allocateRGrid;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::basis;
      using Base::rgrid;
      using Base::phiTot;
      using Base::isAllocatedBasis;
      using Base::isAllocatedRGrid;
      using Base::hasData;
      using Base::isSymmetric;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nBasis;
      using Base::fieldIo;

      /**
      * Calculate the average value of the rgrid_ member.
      */
      double rGridAverage() const override;

   };

   #ifndef RPG_MASK_TPP
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;
   #endif

} // namespace Rpg

#ifndef RPG_MASK_TPP
namespace Prdc {
   extern template class MaskReal< 1, Rpg::FieldIo<1>, Cuda::RField<1> >;
   extern template class MaskReal< 2, Rpg::FieldIo<2>, Cuda::RField<2> >;
   extern template class MaskReal< 3, Rpg::FieldIo<3>, Cuda::RField<3> >;
} 
#endif

} // namespace Pscf
#endif
