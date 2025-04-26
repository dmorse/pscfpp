#ifndef RPC_MASK_H
#define RPC_MASK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"              // parent class template parameter
#include <prdc/cpu/RField.h>      // parent class template parameter
#include <prdc/field/MaskTmpl.h>  // parent class

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * A field to which the total density is constrained.
   *
   * Please refer to the documentation of the base class Prdc::MaskTmpl
   * for more complete API documentation for this class template.
   * The public interface of Rpg::Mask is identical to that of the
   * base class template Prdc::MaskTmpl. 
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Mask : public Prdc::MaskTmpl< D, FieldIo<D>, Prdc::Cpu::RField<D> >
   {

   public:

      /// Base class typedef
      typedef Prdc::MaskTmpl< D, FieldIo<D>, Prdc::Cpu::RField<D> > Base;

      /**
      * Constructor.
      */
      Mask();

      /**
      * Destructor.
      */
      ~Mask();

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

      /**
      * Calculate the average value of the rgrid_ member.
      */
      double rGridAverage() const;

   };

   #ifndef RPC_MASK_TPP
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;
   #endif

} // namespace Rpc

#ifndef RPC_MASK_TPP
namespace Prdc {
   extern template class MaskTmpl< 1, Rpc::FieldIo<1>, Cpu::RField<1> >;
   extern template class MaskTmpl< 2, Rpc::FieldIo<2>, Cpu::RField<2> >;
   extern template class MaskTmpl< 3, Rpc::FieldIo<3>, Cpu::RField<3> >;
} 
#endif

} // namespace Pscf
#endif
