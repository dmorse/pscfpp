#ifndef RPG_W_FIELD_CONTAINER_H
#define RPG_W_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/WFieldsReal.h>     // base class template
#include <prdc/cpu/RField.h>            // template parameter
#include <rpg/field/FieldIo.h>          // template parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * Almost all of the implementation of this class is defined by the base
   * class template Prdc::WFieldsReal . See documentation of that
   * template for information.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class WFieldContainer 
     : public WFieldsReal<D, Prdc::Cuda::RField<D>, Rpg::FieldIo<D> >
   {

   public:

      /// Alias for base class template instantiation
      typedef WFieldsReal<D, RField<D>, FieldIo<D> > Base;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::setNMonomer;
      using Base::allocateRGrid;
      using Base::allocateBasis;
      using Base::allocate;
      using Base::setBasis;
      using Base::setRGrid;
      using Base::readBasis;
      using Base::readRGrid;
      using Base::symmetrize;
      using Base::clear;
      using Base::writeBasis;
      using Base::writeRGrid;
      using Base::basis;
      using Base::rgrid;
      using Base::isAllocatedRGrid;
      using Base::isAllocatedBasis;
      using Base::hasData;
      using Base::isSymmetric;

      /**
      * Set new w fields, in unfolded real-space (r-grid) format.
      *
      * The array fields is an unfolded array that contains fields for
      * all monomer types, with the field for monomer 0 first, etc.
      *
      * \param fields  unfolded array of new w (chemical potential) fields
      */
      void setRGrid(DeviceArray<cudaReal>& fields);

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nBasis;
      using Base::nMonomer;
      using Base::fieldIo;

   private:

      /**
      *  Assign one RField<D> to another: lhs = rhs.
      *
      *  \left lhs  left-hand side of assignment
      *  \left rhs  right-hand side of assignment
      */
      void assignRField(RField<D>& lhs, RField<D> const & rhs) const 
      override;

   };

   #ifndef RPG_W_FIELD_CONTAINER_TPP
   // Suppress implicit instantiation
   extern template class WFieldContainer<1>;
   extern template class WFieldContainer<2>;
   extern template class WFieldContainer<3>;
   #endif

} // namespace Rpg

#ifndef RPG_W_FIELD_CONTAINER_TPP
namespace Prdc {
   // Suppress implicit instantiation
   extern template class WFieldsReal<1, RField<1>, Rpg::FieldIo<1> >;
   extern template class WFieldsReal<2, RField<2>, Rpg::FieldIo<2> >;
   extern template class WFieldsReal<3, RField<3>, Rpg::FieldIo<3> >;
} // namespace Prdc
#endif

} // namespace Pscf
#endif
