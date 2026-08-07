/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/DeviceArray.h>

#include <rp/field/WFieldsBase.tpp>   // base class implementation
#include <rp/field/WFields_u.h>       // class specialziation

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   // Public member function

   /*
   * Set new w-field values, using unfolded array of r-grid fields.
   */
   template <int D>
   void WFields<D,CUT>::setRGrid(DeviceArray<cudaReal>& fields)
   {
      // Create DArray tmp with RField<D,CUT> elements
      DArray< RField<D,CUT> > tmp;
      const int nMonomer = Base::nMonomer();
      tmp.allocate(nMonomer);

      // Associate each RField<D,CUT> with a slice of the unfolded array
      IntVec<D> const & meshDimensions = Base::meshDimensions();
      const int meshSize = Base::meshSize();
      for (int i = 0; i < nMonomer; i++) {
         tmp[i].associate(fields, i*meshSize, meshDimensions);
      }

      // Use tmp array to set w-fields for all monomer types
      bool isSymmetric = false;
      Base::setRGrid(tmp, isSymmetric);
   }

   // Explicit instantiation definitions - base class
   template class WFieldsBase<1,CUT>;
   template class WFieldsBase<2,CUT>;
   template class WFieldsBase<3,CUT>;

   // Explicit instantiation definitions - this class
   template class WFields<1,CUT>;
   template class WFields<2,CUT>;
   template class WFields<3,CUT>;

} // namespace Rp
} // namespace Pscf

