#ifndef RPG_C_FIELD_CONTAINER_H
#define RPG_C_FIELD_CONTAINER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/CFieldsReal.h>   // base class template
#include <rpg/field/FieldIo.h>        // base class template parameter
#include <prdc/cuda/RField.h>         // base class template parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * A list of c fields stored in both basis and r-grid format.
   *
   * This class is simply a named partial specialization of the base 
   * class template Pscf::Prdc::CFieldsReal, designed for use with a GPU.
   * Please see documentation of the base class for API documentation.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class CFieldContainer : public CFieldsReal<D, RField<D>, FieldIo<D> >
   {

   public:

      /// Alias for base class
      typedef CFieldsReal<D, RField<D>, FieldIo<D> >  Base;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::setWriteUnitCell;
      using Base::setNMonomer;
      using Base::allocateRGrid;
      using Base::allocateBasis;
      using Base::allocate;
      using Base::basis;
      using Base::rgrid;
      using Base::writeBasis;
      using Base::writeRGrid;
      using Base::isAllocatedRGrid;
      using Base::isAllocatedBasis;
      using Base::hasData;
      using Base::isSymmetric;
      using Base::setHasData;
      using Base::setIsSymmetric;

   protected:

      using Base::fieldIo;

   };

   // Suppress implicit instantiation
   extern template class CFieldContainer<1>;
   extern template class CFieldContainer<2>;
   extern template class CFieldContainer<3>;

} // namespace Rpg
namespace Prdc {
   // Suppress implicit instantiation of base class
   extern template class CFieldsReal<1, Cuda::RField<1>, Rpg::FieldIo<1> >;
   extern template class CFieldsReal<2, Cuda::RField<2>, Rpg::FieldIo<2> >;
   extern template class CFieldsReal<3, Cuda::RField<3>, Rpg::FieldIo<3> >;
} 
} // namespace Pscf
#endif
