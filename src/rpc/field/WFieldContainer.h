#ifndef RPC_W_FIELD_CONTAINER_H
#define RPC_W_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/WContainerReal.h>     // base class template
#include <prdc/cpu/RField.h>               // template parameter
#include <rpc/field/FieldIo.h>             // template parameter

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * The public interface of this class is identical to that of the base
   * class template Pscf::Prdc::WContainerReal. Please see documentation
   * of that base class for API documentation.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class WFieldContainer 
     : public WContainerReal<D, Prdc::Cpu::RField<D>, Rpc::FieldIo<D> >
   {
   public:

      /// Alias for base class
      typedef WContainerReal<D, RField<D>, FieldIo<D> >  Base;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::setNMonomer;
      using Base::allocateRGrid;
      using Base::deallocateRGrid;
      using Base::allocateBasis;
      using Base::deallocateBasis;
      using Base::allocate;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::symmetrize;
      using Base::clear;
      using Base::basis;
      using Base::rgrid;
      using Base::isAllocatedRGrid;
      using Base::isAllocatedBasis;
      using Base::hasData;
      using Base::isSymmetric;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nBasis;
      using Base::nMonomer;
      using Base::fieldIo;

   private:

      /**
      * Assign one RField<D> to another: lhs = rhs.
      *
      * \param lhs  left-hand side of assignment
      * \param rhs  right-hand side of assignment
      */
      void assignRField(RField<D>& lhs, RField<D> const & rhs) const 
      override;

   };

   #ifndef RPC_W_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class WFieldContainer<1>;
   extern template class WFieldContainer<2>;
   extern template class WFieldContainer<3>;
   #endif

} // namespace Rpc

#ifndef RPC_W_FIELD_CONTAINER_TPP
namespace Prdc {
   // Suppress implicit instantiation
   extern template class WContainerReal<1, RField<1>, Rpc::FieldIo<1> >;
   extern template class WContainerReal<2, RField<2>, Rpc::FieldIo<2> >;
   extern template class WContainerReal<3, RField<3>, Rpc::FieldIo<3> >;
} // namespace Rpc
#endif

} // namespace Pscf
#endif
